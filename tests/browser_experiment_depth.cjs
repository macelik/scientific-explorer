const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:true,executablePath:process.env.CHROMIUM_PATH});
 const page=await browser.newPage({viewport:{width:1500,height:1000}}),errors=[];
 const out=path.resolve(__dirname,'../validation/reload-recursion');
 fs.mkdirSync(out,{recursive:true});
 page.on('pageerror',e=>errors.push(String(e)));
 await page.goto('http://127.0.0.1:8766');
 await page.getByText('300 matching cells',{exact:true}).waitFor({timeout:60000});
 await page.getByRole('button',{name:'TTAGGCTAGGCCGGAA-1 / chr9',exact:true}).click();
 await page.getByRole('button',{name:'X experiments',exact:true}).click();
 const depth=page.getByLabel('Experiment recursion limit',{exact:true});
 await depth.waitFor();assert.equal(await depth.inputValue(),'2');
 await depth.selectOption('4');
 await page.getByRole('button',{name:'Run experiment',exact:true}).click();
 await page.getByRole('heading',{name:'Original → experiment',exact:true}).waitFor({timeout:120000});
 const dl=page.waitForEvent('download');
 await page.getByRole('button',{name:'Export experiment JSON',exact:true}).click();
 await(await dl).saveAs(path.join(out,'depth4-result.json'));
 const result=JSON.parse(fs.readFileSync(path.join(out,'depth4-result.json'),'utf8'));
 assert.equal(result.request.k,4);assert.equal(result.parameters.k,4);
 assert(result.edited.nodes.some(n=>n.depth>1),'Expected deeper nodes');
 assert(result.edited.nodes.every(n=>n.depth<4));
 await page.getByText('Both children and exact split tests',{exact:true}).click();
 const inspect=page.getByLabel('Child inspection depth',{exact:true});
 await inspect.waitFor();
 const inspectPlots=async()=>{
  const plots=await page.locator('.experiment-children [data-node-path]').evaluateAll(els=>els.map(el=>({path:el.dataset.nodePath,shapes:el.querySelector('.js-plotly-plot').layout.shapes})));
  for(const plot of plots){
   const pairs=['original','edited'].map(key=>({key,tree:result[key],node:result[key].nodes.find(n=>n.side===plot.path)})).filter(p=>p.node);
   const intersects=pairs.some(p=>Math.max(p.node.start,result.intervention.target_start)<Math.min(p.node.end,result.intervention.target_end));
   assert.equal(plot.shapes.some(s=>s.type==='rect'),intersects,`Edit overlap ${plot.path}`);
   for(const p of pairs){
    const label=p.key==='original'?'Original':'Edited';
    const expected=p.tree.boundaries.filter(b=>p.node.start<b&&b<p.node.end).map(b=>result.start_bp[b]/1e6);
    const actual=plot.shapes.filter(s=>s.name?.startsWith(label+' accepted')).map(s=>s.x0);
    assert.deepEqual(actual,expected,`Accepted markers ${label} ${plot.path}`);
   }
   for(const s of plot.shapes.filter(s=>s.type==='rect')){
    assert(s.x0>=result.start_bp[result.intervention.target_start]/1e6);
    assert(s.x1<=result.end_bp[result.intervention.target_end-1]/1e6);
    assert.equal(s.yref,'paper');
   }
  }
  return plots;
 };
 const immediate=await inspectPlots();
 const children=page.locator('.experiment-children');
 // Hide sticky header only during element captures so it does not cover plots.
 await page.addStyleTag({content:'header { position: static !important; }'});
 await children.screenshot({path:path.join(out,'children-depth1.png')});
 await inspect.selectOption('2');const deeper=await inspectPlots();
 await children.screenshot({path:path.join(out,'children-depth2.png')});
 await depth.selectOption('3');
 await page.getByText(/Settings changed\. These are the previous run/).waitFor();
 const sessionDownload=page.waitForEvent('download');
 await page.getByRole('button',{name:'download JSON',exact:true}).click();
 const sessionPath=path.join(out,'depth-session.json');await(await sessionDownload).saveAs(sessionPath);
 await page.reload();await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.locator('input[type=file]').setInputFiles(sessionPath);
 await page.getByRole('heading',{name:'X-value experiments',exact:true}).waitFor();
 assert.equal(await depth.inputValue(),'3');
 await page.getByText(/Settings changed\. These are the previous run/).waitFor();
 assert.deepEqual(errors,[]);
 fs.writeFileSync(path.join(out,'browser-results.json'),JSON.stringify({k:4,original_nodes:result.original.nodes.length,edited_nodes:result.edited.nodes.length,immediate,deeper,stale_depth:true,restored_depth:true,errors},null,2));
 console.log('PASS k=4 run, clipped child highlights, accepted breakpoint coordinates, deeper inspection, stale state and session restore');
 await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
