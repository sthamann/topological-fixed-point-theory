import json
from fractions import Fraction as F
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import search


class CharacterSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=search.calculate()
        cls.rs=search.roots()

    def test_full_E8_root_domain(self):
        self.assertEqual(len(self.rs),240)
        self.assertEqual(sum(all(x%2==0 for x in a) for a in self.rs),112)
        self.assertEqual(sum(all(x%2 for x in a) for a in self.rs),128)
        self.assertEqual(len(set(self.rs)),240)

    def test_blind_source_frame_reads_no_dictionary_or_target_constants(self):
        with patch.object(Path,'read_text',side_effect=AssertionError('Blind classification read a file')), \
             patch.object(search,'D',None),patch.object(search,'R',None),patch.object(search,'LAMBDA',None):
            blind=search.blind_source_classification()
        self.assertFalse(blind['reads_M_or_target_projector'])
        self.assertEqual(blind['lambda1_source'],[str(j-4) for j in range(8)])
        self.assertEqual(blind['case_count'],16)
        self.assertEqual(blind['primitive_order_four_count'],12)
        self.assertEqual(blind['D5_A3_candidate_count'],4)

    def test_blind_source_frame_projectors_agree_only_after_transport(self):
        blind=self.result['blind_original_source_frame_control']
        self.assertEqual({(c['u'],c['v']) for c in blind['candidates']},{(1,0),(1,2),(3,0),(3,2)})
        self.assertTrue(all(c['projector_rank']==5 for c in blind['candidates']))
        comparison=blind['post_classification_M_comparison']
        self.assertTrue(all(c['transported_source_projector_equals_target_moment'] for c in comparison))
        self.assertEqual(sum(c['equals_standard_coordinate_P5'] for c in comparison),2)
        self.assertEqual({(c['m'],c['n']) for c in comparison},{(2,3),(2,7),(6,1),(6,5)})

    def test_group_is_complete_and_order_four_subgroup_matches(self):
        signatures={}
        for m in range(16):
            for n in range(8):
                h=search.add(search.scale(m,search.D),search.scale(F(n,8),search.R))
                signatures[search.phases(h,self.rs)]=search.character_order(h,self.rs)
        self.assertEqual(len(signatures),64)
        quarter={search.phases(tuple(F(x) for x in row['h']),self.rs)
                 for row in self.result['all_order_dividing_four_cases']}
        self.assertEqual(quarter,{key for key,order in signatures.items() if order<=4})
        self.assertEqual(len(quarter),16)
        self.assertEqual(self.result['primitive_order_four_count'],12)

    def test_four_candidates_use_actual_components_and_Dynkin_graphs(self):
        actual={(c['m'],c['n']) for c in self.result['candidates']}
        self.assertEqual(actual,{(2,3),(2,7),(6,1),(6,5)})
        for candidate in self.result['candidates']:
            cs=candidate['components']
            self.assertEqual([(c['type'],c['rank'],c['root_count']) for c in cs],[('A3',3,12),('D5',5,40)])
            for c,degrees in zip(cs,([1,1,2],[1,1,1,2,3])):
                cartan=c['cartan_matrix']
                self.assertEqual(sorted(sum(x==-1 for x in row) for row in cartan),degrees)
                self.assertEqual(search.rank(c['simple_roots_doubled']),c['rank'])
            self.assertEqual(candidate['root_phase_counts'],{'0':52,'1/4':64,'1/2':60,'3/4':64})

    def test_injected_wrong_rank_five_projector_rejected(self):
        c=self.result['candidates'][0]
        h=tuple(F(x) for x in c['h'])
        neutral=[a for a in self.rs if search.phase(h,a)==0]
        correct=[[F(x) for x in row] for row in c['projector_from_neutral_moment']]
        self.assertTrue(search.verify_projector(neutral,correct))
        bad=[[F(int(i==j and i in (0,1,2,3,5))) for j in range(8)] for i in range(8)]
        self.assertEqual(search.rank(bad),5)
        self.assertEqual(search.mm(bad,bad),bad)
        self.assertNotEqual(bad,correct)
        with self.assertRaisesRegex(ValueError,'neutral-root moment'):
            search.verify_projector(neutral,bad)

    def test_wrong_Q_over_four_character(self):
        row=self.result['wrong_character_Q_over_4']
        self.assertEqual([row['root_phase_counts'].get(k,0) for k in ('0','1/4','1/2','3/4')],[128,0,112,0])
        self.assertEqual([(c['type'],c['rank'],c['root_count']) for c in row['components']],[('A1',1,2),('E7',7,126)])

    def test_actual_original_fields_and_bilinear_components(self):
        for c in self.result['candidates']:
            self.assertEqual(c['original_complex_source_field_phase_counts'],{'0':2,'1/4':2,'1/2':2,'3/4':2})
            self.assertEqual(c['original_real_source_field_phase_counts'],{'0':4,'1/4':4,'1/2':4,'3/4':4})
            original=c['original_source_bilinear_neutral']
            self.assertEqual(original['root_count'],20)
            self.assertEqual(original['semisimple_rank'],7)
            self.assertEqual(original['central_Cartan_rank'],1)
            self.assertEqual([(x['type'],x['root_count']) for x in original['actual_connected_components']],
                             [('A1',2)]*4+[('A3',12)])

    def test_periodic_control_keeps_actual_quarter_source_distinct(self):
        cases=self.result['state_and_time_scope']['cases']
        self.assertEqual([c['physical_g2_order'] for c in cases],[4,8,8])
        self.assertEqual([c['alpha'] for c in cases],['0','1/4','3/4'])
        self.assertTrue(all(c['constructed_internal_t']==cases[0]['constructed_internal_t'] for c in cases))
        self.assertNotEqual(cases[0]['physical_g2_root_phase_counts'],cases[1]['physical_g2_root_phase_counts'])
        for c in self.result['candidates']:
            self.assertEqual(c['known_ground_pair_marker_phases'],['0','0'])
        self.assertFalse(self.result['new_H_or_state_chosen'])
        self.assertFalse(self.result['physical_source_selected'])

    def test_replay_and_source_change_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)
            fold=directory/'fold-results.json'
            shutil.copyfile(search.FOLD/'fold-results.json',fold)
            output=directory/'result.json'
            search.search(output,fold,atlas_packet=False)
            self.assertTrue(search.verify(output)['verified'])
            replay=subprocess.run([sys.executable,str(Path(search.__file__).resolve()),'verify','--output',str(output)],
                                  text=True,capture_output=True,check=False)
            self.assertEqual(replay.returncode,0,replay.stderr)
            self.assertTrue(json.loads(replay.stdout)['verified'])
            data=json.loads(fold.read_text());data['explicit_Weyl_matrix'][0][0]='0'
            fold.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'Source changed'):
                search.verify(output)


if __name__=='__main__': unittest.main()
