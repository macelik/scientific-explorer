# Artifact checks

Run from repository root: `python3 /tmp/validate_low_coverage_research.py`. Restore the exact source below to that scratch path if needed. These checks validate saved structure, references and selected report/table consistency, not the scientific methods. Combined inventories preserve agent records verbatim and qualify IDs by task. Shared output hashes are checked; original source files are checked by size and mtime against the manifest (their stored hashes are not recomputed here).

```python
from pathlib import Path
import csv, hashlib, json, re

base = Path('research/2026-10-29-low-coverage')
required = {'name','assay','inputs','outputs','mechanism','reference_requirements','low_coverage_handling','homozygous_deletion','evidence','implementation','fit','limitations','source_ids','uncertain'}
implementation = {'code_url','version_or_commit','license','license_url','compute_notes'}
source_fields = {'id','title','url','kind','accessed_at','locator','supports'}
report = {'date':'2026-10-29','scope':'Artifact structure, current completion state, links and shared export integrity; not scientific validation','tasks':[],'errors':[],'warnings':[]}
combined = json.loads((base/'synthesis/methods_inventory.json').read_text())['records']
combined_sources = json.loads((base/'synthesis/sources_inventory.json').read_text())['sources']
all_methods, all_sources = [], []
for d in sorted(base.glob('[a-f]-*')):
    errors = []
    for name in ['TASK.md','findings.md','methods.json','sources.json','pilot.md','status.json']:
        if not (d/name).is_file(): errors.append('Missing '+name)
    state = json.loads((d/'status.json').read_text())
    methods = json.loads((d/'methods.json').read_text())
    sources = json.loads((d/'sources.json').read_text())
    task = state['task_id']
    if state['status'] != 'complete': errors.append('Task incomplete')
    if state['model'] != 'gpt-6.1-terra': errors.append('Model metadata mismatch')
    if len((d/'findings.md').read_text().split()) < 100: errors.append('Findings unfinished')
    ids = {s['id'] for s in sources}
    if len(ids) != len(sources): errors.append('Duplicate source IDs')
    for i,m in enumerate(methods['methods'],1):
        if required-set(m): errors.append(f'Method{i}: missing fields {sorted(required-set(m))}')
        if implementation-set(m['implementation']): errors.append(f'Method{i}: implementation fields missing')
        if set(m['source_ids'])-ids: errors.append(f'Method{i}: unknown source IDs')
        all_methods.append({'record_id':f'{task}:{i}','task_id':task,'task_directory':d.name,'record':m})
    for s in sources:
        if source_fields-set(s): errors.append('Source fields missing '+s.get('id',''))
        all_sources.append({'source_key':f'{task}:{s["id"]}','task_id':task,'record':s})
    report['tasks'].append({'directory':d.name,'status':state['status'],'method_records':len(methods['methods']),'sources':len(sources),'errors':errors})
    report['errors'].extend([d.name+': '+e for e in errors])
if all_methods != combined: report['errors'].append('Combined method inventory differs from agent records')
if all_sources != combined_sources: report['errors'].append('Combined source inventory differs from agent records')
manifest = json.loads((base/'shared/manifest.json').read_text())
for a in manifest['artifacts']:
    p=Path(a['path'])
    if hashlib.sha256(p.read_bytes()).hexdigest()!=a['sha256']: report['errors'].append('Shared artifact hash mismatch '+p.name)
report['shared_hashes_verified'] = len(manifest['artifacts'])
for a in manifest['inputs']:
    st=Path(a['path']).stat()
    if st.st_size!=a['size'] or st.st_mtime_ns!=a['mtime_ns']:report['errors'].append('Source metadata changed '+a['path'])
report['source_metadata_checked'] = len(manifest['inputs'])
links_checked = 0
for p in base.rglob('*.md'):
    content = re.sub(r'```.*?```','',p.read_text(),flags=re.S)
    for link in re.findall(r'\]\(([^)]+)\)',content):
        if link.startswith(('http://','https://','#','mailto:')):continue
        path=link.split('#',1)[0]
        path=re.sub(r':\d+$','',path)
        if path.startswith('file://'):path=path[7:]
        if not (p.parent/path).exists():report['errors'].append(f'Broken local link {p}: {link}')
        links_checked += 1
report['local_links_checked']=links_checked
with (base/'e-per-cell-evidence/artifacts/zero_support_by_stratum.tsv').open() as f:
    pilot={r['region']:r for r in csv.DictReader(f,delimiter='\t') if r['stratum']=='all'}
expected={'chr7_8.3_10.7':(24,7,217,0),'chr7_51.4_53.5':(19,22,255,0),'chr7_67.4_70.0':(26,10,248,0),'chr7_70.0_73.3':(33,0,144,1),'chr9_20_27':(52,0,116,1),'chr9_22.5_26.6':(25,0,204,0)}
for name,values in expected.items():
    r=pilot[name]
    observed=tuple(int(float(r[k])) for k in ['n_bins','raw_whole_zero_cells','at_whole_zero_cells','at_sum_median'])
    if observed!=values:report['errors'].append('Report pilot numbers mismatch '+name)
    if int(r['cells_at_ge10'])!=0:report['errors'].append('Report AT10 claim mismatch '+name)
report['pilot_regions_cross_checked']=len(expected)
report['json_files_parsed']=0
for p in base.rglob('*.json'):
    json.loads(p.read_text());report['json_files_parsed']+=1
report['method_records']=len(all_methods)
report['source_records']=len(all_sources)
report['unique_source_urls']=len({s['record']['url'].rstrip('/') for s in all_sources})
report['scientific_validation_performed']=False
(base/'synthesis/artifact_validation.json').write_text(json.dumps(report,indent=2)+'\n')
src=Path('/tmp/validate_low_coverage_research.py').read_text()
(base/'synthesis/ARTIFACT_CHECKS.md').write_text('# Artifact checks\n\nRun from repository root: `python3 /tmp/validate_low_coverage_research.py`. Restore the exact source below to that scratch path if needed. These checks validate saved structure, references and selected report/table consistency, not the scientific methods. Combined inventories preserve agent records verbatim and qualify IDs by task. Shared output hashes are checked; original source files are checked by size and mtime against the manifest (their stored hashes are not recomputed here).\n\n```python\n'+src+'\n```\n')
print(json.dumps(report,indent=2))
raise SystemExit(bool(report['errors']))

```
