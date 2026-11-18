"""Run exact source and independently implemented finite checks."""
from tfpt_local.independent_check import main as independent
from tfpt_local.mapping_class import main as mapping_class
from tfpt_local.active_search import restore
from tfpt_local.engine import save_json
independent()
mapping_class()
result=restore('results/active-search-final/active-20261024.json')
save_json('results/restore.json',result)
print(result)
