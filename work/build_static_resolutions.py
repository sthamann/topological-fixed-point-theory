import ast, collections, hashlib, json, sys, textwrap
from pathlib import Path
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(repo/'semantic_explorer'))
import extract
corp=extract.Corpus(repo)
baseline=json.loads(Path('/Users/stefanhamann/Documents/Codex/2026-10-13/wie-x20/work/unknown89-baseline.json').read_text())
records=[]
for group in baseline['groups']:
    path=repo/group['script_path']
    text,digest=corp.read(path)
    lines=text.splitlines()
    tree=ast.parse(text)
    script_id='script:'+path.stem
    script_node=corp.script_nodes.get(script_id)
    if not script_node:
        row=corp.by_id[group['claim_id']]
        corp.packet(row)
        script_node=corp.script_nodes[script_id]
    imports=script_node['attrs']['static_inputs'].get('imports',[])
    functions={n.name:n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    for c in group['candidates']:
        line=c['line']
        scopes=[]
        for n in ast.walk(tree):
            if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)) and n.lineno<=line<=n.end_lineno:
                scopes.append((n.end_lineno-n.lineno,n.name,n.lineno,n.end_lineno))
        scope=min(scopes)[1:] if scopes else ('<module>',1,len(lines))
        try: expr=ast.parse(textwrap.dedent(c['value']).strip(),mode='eval').body
        except SyntaxError: expr=None
        names=sorted({n.id for n in ast.walk(expr) if isinstance(n,ast.Name)}) if expr else []
        calls=sorted({ast.unparse(n.func) for n in ast.walk(expr) if isinstance(n,ast.Call)}) if expr else []
        calls_local=[x for x in calls if x.split('.')[-1] in functions]
        calls_imported=[]
        import_names={name for i in imports for name in i.get('names',[])}
        calls_imported=[x for x in calls if x.split('.')[-1] in import_names]
        literal_only=bool(expr) and all(isinstance(n,(ast.Constant,ast.Tuple,ast.List,ast.Set,ast.Dict,ast.UnaryOp,ast.USub,ast.UAdd,ast.Load)) for n in ast.walk(expr))
        if c['kind']=='function_default':
            origin='function_default_in_source'
            origin_note='The expression is the literal default declared in this function signature; runtime overrides are possible.'
        elif '__file__' in c['value']:
            origin='self_source_file_read'
            origin_note='The expression reads the current Python source file (__file__); it is a source-integrity/readback input, not an independent physical datum.'
        elif calls_local:
            origin='helper_expression_in_source'
            origin_note='The expression calls a function defined in this same source file; its output is local code computation, not an external physical input.'
        elif c['kind']=='external_or_state_input':
            if '__file__' in c['value']:
                origin='self_source_file_read'
                origin_note='The expression reads the current Python source file (__file__); it is a source-integrity/readback input, not an independent physical datum.'
            else:
                origin='runtime_or_external_read'
                origin_note='The syntax includes an external/state read; the actual runtime value is not established by source inspection alone.'
        elif literal_only:
            origin='literal_in_source'
            origin_note='The value is explicitly written as a literal/collection in this source file; this establishes code origin, not a derivation from TFPT axioms or physical data.'
        elif calls_local or calls_imported:
            origin='helper_expression_in_source'
            origin_note='The assignment invokes a named local/imported helper visible in the source inventory; source code defines the computation path, but no external physical derivation is implied.'
        else:
            origin='expression_in_source'
            origin_note='The right-hand side is a source-level expression. Its names/functions identify dependencies, but the mapping does not claim that the expression is forced by TFPT or physical evidence.'
        context_start=max(1,line-2);context_end=min(len(lines),line+2)
        line_ref=corp.source(path,line,line)
        context_ref=corp.source(path,context_start,context_end)
        call_chain_sources=[]
        for function_name in calls_local:
            fn=functions[function_name]
            helper_ref=corp.source(path,fn.lineno,fn.end_lineno)
            call_chain_sources.append({'function':function_name,'path':group['script_path'],'line_start':fn.lineno,'line_end':fn.end_lineno,
                                       'sha256':digest,'source_ref':helper_ref['id']})
        try:
            target_node=ast.parse(c['name'],mode='eval').body
            def target_root(node):
                if isinstance(node,(ast.Subscript,ast.Attribute)):return target_root(node.value)
                if isinstance(node,(ast.Tuple,ast.List,ast.Set)):return ','.join(target_root(x) for x in node.elts)
                if isinstance(node,ast.Name):return node.id
                return ast.unparse(node)
            target_roots=[target_root(target_node)]
        except SyntaxError:
            target_roots=[c['name']]
        target_roots=target_roots or names
        uses=[]
        for node in ast.walk(tree):
            if not hasattr(node,'lineno') or node.lineno<=line or not isinstance(node,ast.Name) or not isinstance(node.ctx,ast.Load):
                continue
            if node.id in target_roots:
                nline=node.lineno;txt=lines[nline-1].strip()
                if not uses or uses[-1]['line']!=nline:uses.append({'line':nline,'text':txt})
                if len(uses)>=20:break
        packet=corp.packet(corp.by_id[group['claim_id']])
        claim_sources=[{'id':src['id'],'path':src['path'],'sha256':src['sha256'],'line_start':src['line_start'],'line_end':src['line_end'],'excerpt':src.get('text','')}
                       for src in packet['sources'] if src.get('path','').endswith(('.tex','.md','.txt','.json')) and src.get('path')!='verification/status_ledger.csv']
        inline_comment=lines[line-1].split('#',1)[1].strip() if '#' in lines[line-1] else ''
        preceding_comments=[]
        for comment_line in range(line-1,max(0,line-4),-1):
            stripped=lines[comment_line-1].strip()
            if stripped.startswith('#'):
                preceding_comments.insert(0,{'line':comment_line,'text':stripped})
            elif not stripped:
                continue
            else:
                break
        simple_target=c['name'].split('.')[-1].split('[')[0].strip('() ')
        if c['kind']=='function_default':
            semantic_role='function_parameter_default'
            purpose=f"Default value for parameter `{c['name']}`; it applies when a caller omits that argument. The default declaration alone does not establish whether the value was used or why it was selected."
        elif '__file__' in c['value']:
            semantic_role='source_integrity_readback'
            purpose='Reads this script’s own source text for a static AST firewall/integrity check; it is not a measured or physical input.'
        elif c['candidate_id'].split(':',1)[1].split('[')[0].split('.')[-1] in {'N_CHK','n_eval','total','cnt','si','sii','worst_n','max_use','max_adiff','n_close','n_edge','n_str','coll','ncyc','s_tot','r_tot','c2_fail','n_ref','dev_max','dev_atom','nd_worst','r1_worst','sl_worst','frac_min','dev_fac','dev_z6','conv_dev','id_max','ib_last_fail','N_CROSS','worst_idw'}:
            semantic_role='counter_or_diagnostic'
            purpose='Initializes or updates a local count/diagnostic accumulated by the enclosing computation; this value records the code path rather than supplying independent evidence.'
        elif calls_local:
            semantic_role='local_helper_result'
            purpose=f"Stores the return value of local helper(s) {', '.join(calls_local)}; exact helper definitions are in `call_chain_sources`."
        elif any(isinstance(n,ast.Compare) for n in ast.walk(expr)) if expr else False:
            semantic_role='computed_check'
            purpose='Computes a Boolean consistency condition from local values; it validates this specific implementation and does not establish a general result.'
        elif simple_target in {'N_CHK','n_eval','total','cnt','si','sii','worst_n','max_use','max_adiff','n_close','n_edge','n_str','coll','ncyc','s_tot','r_tot','c2_fail','n_ref','dev_max','dev_atom','nd_worst','r1_worst','sl_worst','frac_min','dev_fac','dev_z6','conv_dev','id_max','ib_last_fail','N_CROSS','worst_idw'}:
            semantic_role='counter_or_diagnostic'
            purpose='Initializes or updates a local count/diagnostic accumulated by the enclosing computation; this value records the code path rather than supplying independent evidence.'
        elif '[' in c['name'] or '(' in c['name']:
            semantic_role='finite_structure_entry'
            purpose=f"Sets an indexed entry in `{target_roots[0] if target_roots else c['name']}` while constructing the finite object in the enclosing function; the value is part of that local representation, not external evidence."
        elif c['name'] in {'mp.mp.dps','mp.dps'}:
            semantic_role='numerical_working_precision'
            purpose=f"Sets mpmath working precision to `{c['value']}` digits for the immediately following numerical calculation. This controls numerical resolution, not a physical parameter."
        elif c['name'].isupper() or simple_target.isupper():
            semantic_role='explicit_code_parameter_or_constant'
            purpose=f"Explicit code-level setting `{c['name']} = {c['value']}`. Its downstream uses are listed; absent a specific supporting citation below, this source does not derive why the setting should be chosen."
        else:
            semantic_role='local_initializer_or_formula'
            purpose=f"Local initializer/formula `{c['name']}` inside `{scope[0]}`; it feeds the local code path shown in `downstream_uses`, without implying a physical derivation."
        statement=''
        statement_start=line;statement_end=line;operation='assignment'
        for candidate_node in ast.walk(tree):
            if not isinstance(candidate_node,(ast.Assign,ast.AnnAssign,ast.AugAssign)) or candidate_node.lineno!=line:continue
            targets=candidate_node.targets if isinstance(candidate_node,ast.Assign) else [candidate_node.target]
            if any(ast.unparse(t)==c['name'] for t in targets):
                statement=ast.get_source_segment(text,candidate_node) or lines[line-1]
                statement_start=candidate_node.lineno;statement_end=candidate_node.end_lineno
                operation='augmented_assignment' if isinstance(candidate_node,ast.AugAssign) else 'assignment'
                break
        if operation=='augmented_assignment':
            semantic_role='algorithmic_in_place_update'
            purpose=f"Applies the in-place operator to `{c['name']}` as part of `{scope[0]}`: `{statement.strip()}`. The update encodes the algorithm shown in this function; it is not an independently selected physical quantity."
        record={
            'id':f"{group['claim_id']}|{script_id}|{c['candidate_id']}",
            'claim_id':group['claim_id'],'script_id':script_id,'candidate_id':c['candidate_id'],
            'candidate_name':c['name'],'candidate_kind':c['kind'],'expression':c['value'],
            'source_origin':origin,'source_origin_note':origin_note,
            'interpretation_status':'source_origin_resolved_role_scoped',
            'interpretation_limit':'This records where the expression comes from in code. It does not establish a physical meaning, an axiom, independent evidence, or a general theorem.',
            'dependencies':names,'calls':calls,'local_call_targets':calls_local,'imported_call_names':calls_imported,
            'script_imports':imports,'call_chain_sources':call_chain_sources,'downstream_uses':uses,
            'claim_source_contexts':claim_sources,'semantic_role':semantic_role,'purpose':purpose,
            'selection_basis':'explicit_in_source' if semantic_role=='explicit_code_parameter_or_constant' else 'algorithmic_or_structural_definition' if semantic_role=='finite_structure_entry' else 'computed_or_diagnostic',
            'inline_comment':inline_comment,'preceding_comments':preceding_comments,
            'selection_rationale':('Source comment: '+inline_comment) if inline_comment else ('Nearby comment: '+' '.join(x['text'] for x in preceding_comments)) if preceding_comments else ('No deeper selection reason is stated in the source; value is a code-level setting.' if semantic_role in {'explicit_code_parameter_or_constant','numerical_working_precision'} else 'Defined by the local algorithm or diagnostic shown in the cited context.'),
            'enclosing_scope':{'name':scope[0],'line_start':scope[1],'line_end':scope[2]},
            'source':{'path':group['script_path'],'sha256':digest,'line':line,'line_text':lines[line-1],
                      'source_ref':line_ref['id'],'context_ref':context_ref['id'],'context_lines':[context_start,context_end],
                      'statement':statement,'statement_lines':[statement_start,statement_end],'operation':operation},
            'method':'source_reviewed_ast_mapping','review_scope':'exact_assignment_and_local_script_context; no script execution'
        }
        records.append(record)
records.sort(key=lambda x:(x['claim_id'],x['script_id'],x['source']['line'],x['candidate_id']))
# Every frozen entry must match the exact current source and be represented once.
assert len(records)==baseline['unknown_expression_count']==195
assert len({x['id'] for x in records})==len(records)
output={'version':1,'repo':str(repo),'baseline_path':'work/unknown89-baseline.json','baseline_group_count':baseline['unknown_group_count'],
        'case_count':len(records),'claim_count':len({x['claim_id'] for x in records}),'script_count':len({x['script_id'] for x in records}),
        'scope':'Resolves the immediate code origin for all expressions grouped in the 89 static_ast unknown entries. Does not infer physical/theory origin.',
        'cases':records}
out=repo/'semantic_explorer/data/static-candidate-resolutions.json'
out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'groups':baseline['unknown_group_count'],'cases':len(records),'claims':output['claim_count'],'scripts':output['script_count'],'origins':dict(collections.Counter(x['source_origin'] for x in records)),'path':str(out)},ensure_ascii=False))
