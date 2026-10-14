#!/usr/bin/env python3
"""Read-only lexical call analysis; never invokes a compiler or changes the input tree.

SPDX-License-Identifier: BSD-3-Clause
"""
from pathlib import Path
from collections import defaultdict
import re, json, subprocess, sys, gzip, base64, argparse
from html import escape as esc

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent,help='Source checkout; defaults to this repository')
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().with_name('icon-land-call-tree.html'))
parser.add_argument('--repository-url',default='',help='Optional public repository URL for source links')
args=parser.parse_args()
ROOT=args.source.resolve()
DEST=args.output.resolve()
files = [x for x in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
         if x.lower().endswith('.f90') and x.startswith(('src/','support/','externals/jsbach/src/','ragnarok/'))]
source = {f:(ROOT/f).read_text(errors='replace').splitlines() for f in files}
nodes, bodies, imports, module_imports, bindings, generics = {}, {}, {}, defaultdict(list), {}, defaultdict(list)
binding_sites={}
by_name, by_module = defaultdict(list), defaultdict(list)
cpp_disabled={'__no_icon_atmo__','__no_icon_ocean__','__no_icon_waves__','__no_icon_testbed__',
              '__no_icon_comin__','__no_icon_coupling__','__no_jsbach_hd__','__no_aes__','__no_nwp__'}

def no_comment(text):
    result, quote, i = '', None, 0
    while i < len(text):
        c=text[i]
        if quote:
            result+=c
            if c==quote:
                if i+1<len(text) and text[i+1]==quote:
                    result+=text[i+1]; i+=1
                else: quote=None
        elif c in "'\"": quote=c; result+=c
        elif c=='!': break
        else: result+=c
        i+=1
    return result

def statements(lines):
    buf, start, guards = '', 0, []
    for number, raw in enumerate(lines,1):
        stripped=raw.strip()
        if stripped.startswith('#'):
            m=re.match(r'#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)',stripped)
            if m:
                op, expr=m.groups(); expr=expr.strip()
                if op in ('ifdef','ifndef','if'): guards.append(op+' '+expr)
                elif op=='endif' and guards: guards.pop()
                elif op in ('else','elif') and guards: guards[-1]='alternative to '+guards[-1]
            continue
        s=no_comment(raw).strip()
        if not s: continue
        if not buf: start=number
        if s.startswith('&'): s=s[1:].lstrip()
        more=s.endswith('&')
        if more: s=s[:-1].rstrip()
        buf+=' '+s
        if not more:
            # Split executable semicolon statements without breaking literals.
            quote=None; part=''; chunks=[]
            for c in buf.strip():
                if quote:
                    part+=c
                    if c==quote: quote=None
                elif c in "'\"": quote=c;part+=c
                elif c==';': chunks.append(part.strip());part=''
                else:part+=c
            if part.strip():chunks.append(part.strip())
            for part in chunks: yield start, number, part, list(guards)
            buf=''
    if buf.strip():yield start,len(lines),buf.strip(),list(guards)

def area(file):
    if file.startswith('externals/jsbach/src/'):
        part=file.split('/')[3]
        return 'land-'+part
    if file.startswith('ragnarok/'):return 'kokkos'
    return 'host'

def parse_use(s):
    m=re.match(r'use\s*(?:,\s*[^:]+::\s*)?(\w+)\s*(?:,\s*only\s*:\s*(.*))?$',s,re.I)
    if not m:return None
    mod=m[1].lower();only=m[2]
    names={}
    if only:
        for x in only.split(','):
            terms=x.strip().split('=>')
            if re.match(r'^\w+$',terms[0].strip()):names[terms[0].strip().lower()]=terms[-1].strip().lower()
    return {'module':mod,'names':names,'all':only is None}

start_re=re.compile(r'^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*(?:(?:integer|real|logical|character|complex|double\s+precision|type|class)(?:\s*\([^)]*\))?\s+)?(subroutine|function|program)\s+(\w+)\b',re.I)
for file in files:
    mod=''
    stack=[]; typ=None; iface=[]
    for first,last,s,guards in statements(source[file]):
        low=s.lower()
        m=re.match(r'^module\s+(?!procedure\b|function\b|subroutine\b)(\w+)',low)
        if m:mod=m[1];continue
        if re.match(r'^end\s*module\b',low):mod='';continue
        if re.match(r'^(?:abstract\s+)?interface\b',low):
            name=re.sub(r'^(?:abstract\s+)?interface\s*','',low).strip()
            iface.append(name);continue
        if re.match(r'^end\s*interface\b',low):
            if iface:iface.pop()
            continue
        if re.match(r'^end\s*type\b',low):typ=None;continue
        m=re.match(r'^type\s*(?:,.*?::\s*|::\s*)(\w+)\b',low)
        if m and not stack:
            typ=m[1];bindings[(mod,typ)]={};continue
        m=re.match(r'^type\s+(\w+)\s*$',low)
        if m and not stack:typ=m[1];bindings[(mod,typ)]={};continue
        if typ and not stack:
            m=re.match(r'^procedure\b.*?::\s*(\w+)\s*(?:=>\s*(\w+))?',low)
            if m:
                bindings[(mod,typ)][m[1]]=m[2] or m[1]
                binding_sites[(mod,typ,m[1])]={'file':file,'line':first,'code':s,'guards':guards}
            continue
        m=re.match(r'^module\s+procedure\s+(.*)',low)
        if m and iface and iface[-1]:
            generics[(mod,iface[-1])]+=[x.strip() for x in m[1].split(',')]
            continue
        if re.match(r'^end\s*(subroutine|function|program)\b',low):
            if stack:nodes[stack.pop()]['end']=last
            continue
        m=start_re.match(s)
        if m:
            kind,name=m[1].lower(),m[2]
            parent=stack[-1] if stack else None
            ident=f'{mod or "global"}::{name.lower()}@{file}:{first}'
            bm=re.search(r'bind\s*\(\s*c\s*(?:,\s*name\s*=\s*[\'"]([^\'"]+)[\'"])?',s,re.I)
            nodes[ident]={'id':ident,'name':name,'module':mod,'file':file,'line':first,'end':last,
                          'kind':'interface' if iface else kind,'area':area(file),'language':'Fortran',
                          'signature':s,'parent':parent,'guards':guards,'bind':bm[1] or name if bm else None,
                          'snippet':'\n'.join(f'{i+1:>5}  {x}' for i,x in enumerate(source[file][first-1:first+10],first-1))}
            imports[ident]=[];bodies[ident]=[];stack.append(ident)
            by_name[name.lower()].append(ident);by_module[(mod,name.lower())].append(ident)
            continue
        u=parse_use(s)
        if u:
            if stack:imports[stack[-1]].append(u)
            else:module_imports[mod].append(u)
        if stack:bodies[stack[-1]].append((first,last,s,guards))

def lookup_module(module,name,seen=None):
    seen=set() if seen is None else seen
    if (module,name) in seen:return []
    seen=seen|{(module,name)}
    ids=by_module.get((module,name),[])
    defs=[x for x in ids if nodes[x]['kind']!='interface']
    if defs:return defs
    if ids:return ids
    if (module,name) in generics:
        return sum((lookup_module(module,n,seen) for n in generics[(module,name)]),[])
    out=[]
    for u in module_imports[module]:
        if name in u['names']:out+=lookup_module(u['module'],u['names'][name],seen)
        elif u['all']:out+=lookup_module(u['module'],name,seen)
    return list(dict.fromkeys(out))

def context_imports(ident):
    out=[];n=nodes[ident]
    while ident:
        out+=imports.get(ident,[]);ident=nodes[ident].get('parent')
    return out+module_imports[n['module']]

def resolve(ident,name,allow_unique=False):
    n=nodes[ident];parent=ident
    # Contained procedures and their own recursion take precedence.
    local=[x for x in by_name.get(name,[]) if nodes[x]['file']==n['file'] and nodes[x]['module']==n['module'] and
           (nodes[x].get('parent')==ident or x==ident or nodes[x].get('parent')==n.get('parent'))]
    if local:return local,'direct'
    imported_mod=None
    for u in context_imports(ident):
        if name in u['names']:
            imported_mod=u['module'];ids=lookup_module(u['module'],u['names'][name])
            if ids:return ids,'generic' if len(ids)>1 else 'direct'
    ids=lookup_module(n['module'],name)
    if ids:return ids,'generic' if len(ids)>1 else 'direct'
    for u in context_imports(ident):
        if u['all']:
            ids=lookup_module(u['module'],name)
            if ids:return ids,'generic' if len(ids)>1 else 'direct'
    candidates=[x for x in by_name.get(name,[]) if nodes[x]['kind']!='interface']
    if allow_unique and len(candidates)==1:return candidates,'inferred'
    return [],imported_mod or ''

def named(name,module=None):
    ids=by_module.get((module,name.lower()),[]) if module else by_name.get(name.lower(),[])
    defs=[x for x in ids if nodes[x]['kind']!='interface']
    ids=defs or ids
    if len(ids)!=1:raise ValueError(f'Expected unique {module}:{name}, got {len(ids)}')
    return ids[0]

edges=defaultdict(list)
def external(name,module='',dynamic=False):
    ident=('dynamic::' if dynamic else 'external::')+module.lower()+'::'+name.lower()
    if ident not in nodes:
        nodes[ident]={'id':ident,'name':name,'module':module,'file':'','line':0,'end':0,
                      'kind':'dynamic' if dynamic else 'external','area':'unresolved' if dynamic else 'external',
                      'language':'unknown','signature':'','guards':[],
                      'note':'Receiver type is not statically resolved.' if dynamic else 'Definition is outside the inspected source, an intrinsic, or a retained build-dependent branch.'}
    return ident

def edge(fr,to,label,line=0,kind='direct',code='',guards=None,file=None):
    key=(to,label,kind)
    old=next((e for e in edges[fr] if (e['to'],e['label'],e['kind'])==key),None)
    site={'file':file or nodes[fr].get('file',''),'line':line,'code':code[:240],'guards':guards or []}
    if old:
        if len(old['sites'])<2:old['sites'].append(site)
        old['count']+=1
    else:edges[fr].append({'to':to,'label':label,'kind':kind,'sites':[site],'count':1})

# Only these receiver bindings are promoted to resolved targets; others stay dynamic.
method_modules={'run_tasks':'mo_jsb_model_class','set_options':'mo_jsb_model_class','build_task_schedule':'mo_jsb_model_class',
 'queue_task':'mo_jsb_model_class','build_tile_lists':'mo_jsb_model_class','configure_processes':'mo_jsb_model_class',
 'do_it':'mo_jsb_task_class','aggregate_weighted_by_fract_2d':'mo_jsb_tile_class'}
method_names={'configure_processes':'configure_model_processes'}
known_function_names={k for k,v in by_name.items() if any(nodes[x]['kind']=='function' for x in v)}
for ident in list(bodies):
    if nodes[ident]['kind']=='interface':continue
    for first,last,s,guards in bodies[ident]:
        cm=re.search(r'\bcall\s+([\w]+(?:\s*%\s*[\w]+)*)',s,re.I)
        call_name=None
        if cm:
            call_name=re.sub(r'\s+','',cm[1]).lower()
            if '%' in call_name:
                method=call_name.split('%')[-1]
                if method in method_modules:
                    try:ids=[named(method_names.get(method,method),method_modules[method])]
                    except ValueError:ids=[]
                    typ='binding'
                else:ids=[];typ='dynamic'
                if not ids:ids=[external(call_name,nodes[ident]['module'],True)]
            else:
                ids,typ=resolve(ident,call_name,True)
                if not ids:ids=[external(call_name,typ)];typ='external'
            for target in ids:edge(ident,target,cm[1],first,typ,s,guards)
        clean=re.sub(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"",'',s)
        for match in re.finditer(r'(?<![\w%])(\w+)\s*\(',clean):
            name=match[1].lower()
            if name==call_name or name not in known_function_names|{'mm_initialize','mm_finalize'}:continue
            ids,typ=resolve(ident,name)
            if not ids and name in {'mm_initialize','mm_finalize'}:ids=[external(name,'memman')];typ='external'
            for target in ids:edge(ident,target,name,first,'function' if typ=='direct' else typ,s,guards)

# Source-verified virtual dispatch targets, resolved through each task constructor.
tasks=[]; init_id=named('jsbach_init','mo_jsb_model_init')
for first,last,s,guards in bodies[init_id]:
    m=re.search(r'call\s+model%queue_task\s*\(\s*(\w+)\s*,\s*[\'"]([^\'"]+)',s,re.I)
    if not m:continue
    process,task=m.groups();ctor=named('create_task_'+task)
    body=' '.join(x[2] for x in bodies[ctor])
    tm=re.search(r'allocate\s*\(\s*(tsk_\w+)\s*::',body,re.I)
    if not tm:raise ValueError('Cannot resolve constructor '+task)
    typ=tm[1].lower();mod=nodes[ctor]['module'];b=bindings[(mod,typ)]
    tasks.append({'name':task,'process':process.upper(),'order':len(tasks)+1,'queueLine':first,
                  'constructor':ctor,'type':typ,'integrate':named(b['integrate'],mod),
                  'aggregate':named(b['aggregate'],mod),
                  'integrateSite':dict(binding_sites[(mod,typ,'integrate')]),
                  'aggregateSite':dict(binding_sites[(mod,typ,'aggregate')])})
assert len(tasks)==22
doit=named('do_it','mo_jsb_task_class')
for action,verb in [('integrate','Integrate'),('aggregate','Aggregate')]:
    target=external('this%'+verb,'mo_jsb_task_class',True)
    sites=[e for e in edges[doit] if e['to']==target]
    at=sites[0]['sites'][0]['line'] if sites else nodes[doit]['line']
    for t in {t['name']:t for t in tasks}.values():
        site=t[action+'Site']
        edge(target,t[action],t['process']+' '+t['name'],site['line'],'dispatch',site['code'],site['guards'],site['file'])

def cpp(ident,name,file,pattern,note=''):
    lines=(ROOT/file).read_text().splitlines()
    line=next(i+1 for i,s in enumerate(lines) if re.search(pattern,s))
    nodes[ident]={'id':ident,'name':name,'module':name.rsplit('::',1)[0] if '::' in name else 'extern C',
                  'file':file,'line':line,'end':line,'kind':'C++ function','area':'kokkos','language':'C++',
                  'signature':lines[line-1].strip(),'guards':[],'note':note,
                  'snippet':'\n'.join(f'{i+1:>5}  {x}' for i,x in enumerate(lines[line-1:line+18],line-1))}
    return ident

rf='ragnarok/support/icon_bridge/ragnarok.cpp'
cf='ragnarok/support/common/execution_space.hpp'
ci=cpp('cpp::init_c','init_ragnarok_c',rf,r'^void init_ragnarok_c\(')
ri=cpp('cpp::init','ragnarok::init',rf,r'^void init\(')
cc=cpp('cpp::finalize_c','finalize_ragnarok',rf,r'^void finalize_ragnarok\(')
rc=cpp('cpp::finalize','ragnarok::finalize',rf,r'^void finalize\(')
create=cpp('cpp::exec','icon::create_default_execution_space',cf,r'^inline void create_default_execution_space\(std::intptr_t')
destroy=cpp('cpp::exec_destroy','icon::destroy_default_execution_space',cf,r'^inline void destroy_default_execution_space\(')
bridge=cpp('cpp::bridge','icon::init_ragnarok_f2c','ragnarok/support/icon_bridge/icon_bridge.cpp',r'extern "C" void init_ragnarok_f2c')
domain=cpp('cpp::domain','icon::domain::init','ragnarok/support/icon_bridge/icon_domain.cpp',r'^void init\(')
timer=cpp('cpp::timer_c','init_ragnarok_timer','ragnarok/support/icon_bridge/icon_timer.cpp',r'^extern "C" void init_ragnarok_timer')
timer_backend=cpp('cpp::timer','timer::connect_icon_backend','ragnarok/support/icon_bridge/icon_timer.cpp',r'^void connect_icon_backend')
for fname,to in [('init_ragnarok_c',ci),('finalize_ragnarok',cc),('init_ragnarok_f2c',bridge),('init_ragnarok_timer',timer)]:
    candidates=[n for n in nodes.values() if n.get('bind')==fname]
    for n in candidates:edge(n['id'],to,'BIND(C) → extern C',n['line'],'bridge',n['signature'])
for fr,to,label in [(ci,ri,'ragnarok::init'),(cc,rc,'ragnarok::finalize'),(ri,create,'bind to OpenACC stream'),
                    (rc,destroy,'release execution space'),(bridge,domain,'domain::init'),(timer,timer_backend,'connect_icon_backend')]:
    edge(fr,to,label,nodes[fr]['line'],'cpp')
for fr,name in [(ri,'Kokkos::initialize'),(rc,'Kokkos::finalize'),(timer_backend,'Timer::set_backend')]:
    edge(fr,external(name,'C++ library'),name,nodes[fr]['line'],'external')
for fun in ['get_domain_info','get_patch_info']:
    edge(domain,named(fun,'mo_ragnarok_support'),'callback: '+fun,nodes[domain]['line'],'callback',
         'Callback through the function table registered by init_ragnarok_support; may be reached via C++ patch helpers.')

notes={
 'jsbach_model':'Standalone land driver. Sets up MPI work groups, models, grid, tiles, forcing, output and restart; loops over timesteps and model instances; finalizes resources.',
 'run_one_timestep':'For one model: start timestep, prepare forcing, loop over local nproma blocks, set chunk options, run the scheduled tasks, finish timestep. CPU block parallelism is on the Fortran side here.',
 'jsbach_init':'Queues the scientific tasks, initializes tile/process state and constructs the reusable task schedule. The queue contains 22 entries and 21 distinct task names.',
 'build_task_schedule':'Walks each task through the tile tree. SKIP/INHERIT branches are omitted; children are integrated before parent aggregation. Inactive processes are omitted at runtime.',
 'run_tasks':'For each queued active task, executes its precomputed tile/action steps. Optional multiple OpenACC queues are joined before aggregation and at task completion.',
 'do_it':'Dynamic task dispatch: TASK_INTEGRATE calls the concrete Integrate binding; TASK_AGGREGATE calls Aggregate. Shown candidate targets are source-resolved from queued task constructors, not all executed at every call.',
 'init_ragnarok':'Initializes Kokkos on the GPU already selected by OpenACC, sharing its queue-1 stream. The CPU path passes no device or stream.',
 'init_ragnarok_support':'Registers host function pointers for grid, halo exchange, logging and timers, then connects the timer backend.',
 'calc_vol_heat_capacity':'Scalar soil-layer formula; the properties task invokes it over cells and soil levels. A contained first kernel for a future C++ port.',
 'calc_soil_temperature':'Soil-column numerical calculation. Preserve sequential dependencies along the vertical axis when introducing parallel execution.',
}
for n in nodes.values():
    if n['name'].lower() in notes:n['note']=notes[n['name'].lower()]

def tree_node(name=None,module=None,label=None,note='',children=None,ident=None):
    ident=ident or named(name,module)
    return {'id':ident,'label':label or nodes[ident]['name'],'note':note,'children':children or []}

def phase(label,children,note=''):
    return {'label':label,'note':note,'children':children,'phase':True}

lifecycle=tree_node('icon',label='icon · program entry',children=[
 phase('01 · Master setup',[
    tree_node('start_mpi','mo_mpi'),tree_node('init_master_control','mo_master_init'),
    tree_node(ident=external('mm_initialize','memman'),label='mm_initialize() · MemMan')]),
 tree_node('jsbach_model','mo_jsbach_model',label='02 · jsbach_model · standalone driver',children=[
    phase('Initialization',[
      tree_node('read_infrastructure_namelists_for_icon_land','mo_read_namelists_iface'),
      tree_node('set_mpi_work_communicators','mo_mpi'),tree_node('init_ragnarok','mo_ragnarok'),
      tree_node('jsbach_setup_models','mo_jsb_base'),tree_node('build_decomposition','mo_build_decomposition',note='Alternative: build_sitelevel_decomposition for site-level runs.'),
      tree_node('init_ragnarok_support','mo_ragnarok_support'),tree_node('setup_forcing','mo_jsb4_forcing'),
      tree_node('setup_grid_patch','mo_jsb_model_init',label='jsbach_setup_grid → setup_grid_patch'),tree_node('jsbach_setup_tiles','mo_jsb_base'),
      tree_node('jsbach_init','mo_jsb_model_init'),tree_node('init_forcing','mo_jsb4_forcing'),
      phase('Restart restore · conditional',[tree_node('jsbach_init_after_restart','mo_jsb_model_init')]),
      tree_node('init_name_list_output','mo_name_list_output_init')]),
    phase('Every timestep · every model',[
      tree_node('run_one_timestep','mo_jsbach_model',children=[
        tree_node('jsbach_start_timestep','mo_jsb_interface'),tree_node('get_standalone_driver','mo_jsb4_forcing'),
        phase('Every local block · OpenMP on CPU',[
          tree_node('set_options','mo_jsb_model_class',label='model%Set_options'),
          tree_node('run_tasks','mo_jsb_model_class',label='model%Run_tasks',children=[
            tree_node('do_it','mo_jsb_task_class',label='task%Do_it · integrate or aggregate',children=[
              tree_node(ident=external('this%Integrate','mo_jsb_task_class',True),label='this%Integrate · concrete update routines'),
              tree_node(ident=external('this%Aggregate','mo_jsb_task_class',True),label='this%Aggregate · concrete aggregation routines')])])]),
        tree_node('jsbach_finish_timestep','mo_jsb_interface')]),
      tree_node('update_statistics','mo_derived_variable_handling'),tree_node('write_name_list_output','mo_name_list_output'),
      phase('Restart events · conditional',[
        tree_node(ident=external('restartdescriptor%writerestart','mo_jsbach_model',True),label='restartDescriptor%writeRestart',note='Runtime binding on restartDescriptor; implementation is configuration-dependent.')])]),
    phase('Finalization',[
      tree_node('finalize_external_forcing','mo_jsb4_forcing'),tree_node('jsbach_finalize','mo_jsb_model_final'),
      tree_node('finalize_ragnarok','mo_ragnarok'),tree_node('close_name_list_output','mo_name_list_output')])]),
 phase('03 · Program shutdown',[
    tree_node(ident=external('mm_finalize','memman'),label='mm_finalize() · MemMan'),tree_node('stop_mpi','mo_mpi')])])

bridge_tree=phase('Fortran ↔ C++ / Kokkos',[
 tree_node('init_ragnarok','mo_ragnarok',label='Runtime initialization'),
 tree_node('init_ragnarok_support','mo_ragnarok_support',label='Host services / callback table'),
 tree_node('finalize_ragnarok','mo_ragnarok',label='Runtime finalization')],
 'Runtime and host-service interop currently present in the source. No C++ land-physics dispatch is connected yet.')

# Call-site excerpts aid review without shipping the full repository in the HTML.
for ident,n in nodes.items():
    if n.get('file') and n['language']=='Fortran':
        n['end']=max(n['end'],n['line'])
    n.pop('parent',None)
incoming=defaultdict(list)
for fr,es in edges.items():
    for e in es:incoming[e['to']].append({'from':fr,'kind':e['kind'],'label':e['label'],'line':e['sites'][0]['line']})
snapshot=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
payload={'repo':'.','repositoryUrl':args.repository_url.rstrip('/'),'commit':snapshot,'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
         'nodes':nodes,'edges':dict(edges),'incoming':dict(incoming),'tasks':tasks,'lifecycle':lifecycle,'bridge':bridge_tree,
         'stats':{'files':len(files),'routines':sum(n['language']=='Fortran' and n['kind']!='interface' for n in nodes.values()),
                  'edges':sum(len(x) for x in edges.values()),'landFiles':sum(f.startswith('externals/jsbach/src/') for f in files)},
         'methodology':[
           'Static source analysis of tracked Fortran declarations, CALL statements, known function references, USE imports, renames and selected type-bound bindings. No compiler, application build, or runtime tracing.',
           'The application view is a source-verified overview of the retained standalone land path. Phase headings are editorial groups, not additional subroutines. Task ordering comes directly from jsbach_init.',
           'The routine explorer includes raw source branches. Preprocessor guards are recorded, not evaluated. Dormant branches for removed ICON components can therefore appear there.',
           'Generic calls show possible overloads. Dynamic task targets are resolved from constructor types and Integrate/Aggregate bindings; runtime process/tile activation still controls execution.',
           'Other polymorphic calls, procedure-pointer calls, macro-generated calls and unresolved dependencies are incomplete. Unresolved nodes are labeled; uninitialized external submodules are not scanned.',
           'Function references are lexical candidates, not compiler-proven calls. Edges marked inferred use a unique routine-name match rather than complete Fortran name resolution.',
           'The C++ bridge is traced manually for initialization, finalization and host callbacks. Other C++ implementation call graphs are not comprehensively indexed.',
           'Repeated callees unfold as separate tree occurrences. Recursion is stopped at a cycle marker. Static source order is not a measured execution sequence or frequency.'
         ]}
# Compact the self-contained file; hydrate source filenames and reverse edges in the browser.
payload.pop('incoming')
ids={old:str(i+1) for i,old in enumerate(nodes)}
file_table=sorted({n['file'] for n in nodes.values() if n.get('file')} | {s['file'] for es in edges.values() for e in es for s in e['sites'] if s.get('file')})
file_index={f:i for i,f in enumerate(file_table)}
compact_nodes={}
for old,n in nodes.items():
    n['id']=ids[old]
    n['file']=file_index[n['file']] if n.get('file') else -1
    n.pop('signature',None);n.pop('bind',None)
    compact_nodes[ids[old]]=n
compact_edges={}
for old,es in edges.items():
    for e in es:
        e['to']=ids[e['to']]
        for site in e['sites']:site['file']=file_index.get(site['file'],-1)
    compact_edges[ids[old]]=es
def compact_tree(t):
    if t.get('id'):t['id']=ids[t['id']]
    for c in t.get('children',[]):compact_tree(c)
for t in tasks:
    for key in ['constructor','integrate','aggregate']:t[key]=ids[t[key]]
    for key in ['integrateSite','aggregateSite']:t[key]['file']=file_index[t[key]['file']]
compact_tree(lifecycle);compact_tree(bridge_tree)
payload['nodes']=compact_nodes;payload['edges']=compact_edges;payload['files']=file_table
html=Path(__file__).resolve().with_name('icon-land-call-tree-template.html').read_text()
utility=re.compile(r'^(message|finish|timer_|new_timer|debug_on|timer_on|deallocate|datetime|timedelta|get_model$|get_grid$|int2string$|real2string$|assign_if_present$|get_omp_thread$|is_omp_inside_serial$|c_f_pointer$|getptstring)',re.I)
native_ids={i for i,n in compact_nodes.items() if n['area'].startswith('land-')}
for i in list(native_ids):native_ids.update(e['to'] for e in compact_edges.get(i,[]))
def gather(t):
    if t.get('id'):native_ids.add(t['id'])
    for child in t.get('children',[]):gather(child)
gather(lifecycle);gather(bridge_tree)
for i in list(native_ids):
    if compact_nodes[i]['area']=='kokkos':native_ids.update(e['to'] for e in compact_edges.get(i,[]))
def native_link(i,label=None):
    name=esc(label or compact_nodes[i]['name'])
    if i not in native_ids:return '<span>'+name+'</span>'
    return '<a class="native-link" href="#routine-'+i+'">'+name+'</a>'
def native_tree(item,depth=0,seen=None):
    seen=set() if seen is None else seen
    ident=item.get('id');n=compact_nodes.get(ident,{})
    label=item.get('label',n.get('name',''))
    title=native_link(ident,label) if ident else esc(label)
    children=item.get('children',[])
    cycle=ident in seen
    if not children and ident and depth<2 and not cycle:
        children=[{'id':e['to'],'label':e['label']} for e in compact_edges.get(ident,[]) if not utility.match(compact_nodes[e['to']]['name'])][:8]
    note='<div class="node-sub">'+esc(item['note'])+'</div>' if item.get('note') else ''
    if cycle:children=[];note+='<span class="badge dynamic">cycle</span>'
    if children:
        content=''.join(native_tree(c,depth+1,seen|({ident} if ident else set())) for c in children)
        return '<li><details'+(' open' if depth<2 or item.get('phase') else '')+'><summary class="tree-row"><span class="caret">›</span>'+title+'</summary>'+note+'<ul class="tree-list">'+content+'</ul></details></li>'
    return '<li><div class="tree-row">'+title+'</div>'+note+'</li>'
application='<ul class="tree-list">'+native_tree(lifecycle)+'</ul>'
native_tasks=[]
for t in tasks:
    native_tasks.append('<li><details><summary class="tree-row"><span class="caret">›</span><span class="task-num">'+str(t['order']).zfill(2)+'</span><span class="node-text">'+esc(t['name'])+'</span><span class="badge land">'+esc(t['process'])+'</span></summary><ul class="tree-list"><li><div class="tree-row">Integrate → '+native_link(t['integrate'])+'</div></li><li><div class="tree-row">Aggregate → '+native_link(t['aggregate'])+'</div></li></ul></details></li>')
catalog=[]
for ident in sorted(native_ids,key=lambda i:(not compact_nodes[i]['area'].startswith('land-'),compact_nodes[i]['name'].lower())):
    n=compact_nodes[ident];file=file_table[n['file']] if n['file']>=0 else ''
    calls=[e for e in compact_edges.get(ident,[]) if not utility.match(compact_nodes[e['to']]['name'])]
    links=' '.join(native_link(e['to']) for e in calls[:8])
    tail='<div class="native-meta">Further host calls are available in the enhanced explorer.</div>' if len(calls)>8 else ''
    snippet='\n'.join(n.get('snippet','').splitlines()[:4])
    content='<div class="native-meta">'+esc(file+':'+str(n['line']) if file else 'External or unresolved target')+'</div>'
    if n.get('note'):content+='<p class="detail-note">'+esc(n['note'])+'</p>'
    if snippet:content+='<pre>'+esc(snippet)+'</pre>'
    if links:content+='<div class="native-calls">Calls / candidates: '+links+'</div>'+tail
    catalog.append('<details class="source-card" id="routine-'+ident+'"><summary><code>'+esc(n['name'])+'</code> <small>'+esc(n.get('module',''))+'</small></summary>'+content+'</details>')
detail_id=lifecycle['children'][1]['children'][1]['children'][0]['id']
dn=compact_nodes[detail_id]
details='<h2 class="detail-name">'+esc(dn['name'])+'</h2><p class="detail-note">'+esc(dn.get('note',''))+'</p><p>'+native_link(detail_id,'View routine reference')+'</p><p class="fallback-tip">The tree and reference links work without JavaScript. Enhanced search appears when available.</p>'
replacements={
 '__NATIVE_APPLICATION__':application,'__NATIVE_TASKS__':'<ul class="tree-list">'+''.join(native_tasks)+'</ul>',
 '__NATIVE_BRIDGE__':'<ul class="tree-list">'+native_tree(bridge_tree)+'</ul>',
 '__NATIVE_ROUTINES__':''.join(catalog),'__NATIVE_DETAILS__':details,
 '__NATIVE_METHOD__':''.join('<li>'+esc(x)+'</li>' for x in payload['methodology']),
 '__NATIVE_COUNTS__':''.join('<span><strong>'+str(v)+'</strong>'+label+'</span>' for v,label in [(payload['stats']['files'],'Fortran files indexed'),(payload['stats']['routines'],'routine definitions'),(22,'queued tasks')]),
 '__BRANCH__':esc(payload['branch']),'__COMMIT__':snapshot[:7],'__REPO__':esc(args.repository_url or 'Source checkout: .'),
}
for key,value in replacements.items():html=html.replace(key,value)
encoded=base64.b64encode(gzip.compress(json.dumps(payload,separators=(',',':'),ensure_ascii=False).encode('utf-8'),compresslevel=9,mtime=0)).decode('ascii')
DEST.with_name('icon-land-call-tree-data.js').write_text('window.ICON_LAND_CALL_DATA = '+json.dumps(encoded)+';\n')
DEST.write_text(html)
print(json.dumps({'output':str(DEST),'bytes':DEST.stat().st_size,'files':len(files),**payload['stats'],'taskEntries':len(tasks),'distinctTasks':len({t['name'] for t in tasks})},indent=2))
