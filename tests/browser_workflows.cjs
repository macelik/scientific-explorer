const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const browser=await chromium.launch({headless:true,executablePath:process.env.CHROMIUM_PATH});
 const page=await browser.newPage({viewport:{width:1500,height:1080}}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 await page.goto('http://127.0.0.1:8766', {waitUntil:'domcontentloaded',timeout:60000});
 await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'TTAGGCTAGGCCGGAA-1 / chr9',exact:true}).click();
 await page.locator('.explorer .node-card .js-plotly-plot').first().waitFor();
 await page.getByRole('button',{name:'Select window',exact:true}).click();
 const coords=async(selector,x0,x1)=>{
  const plot=page.locator(selector).first();await plot.scrollIntoViewIfNeeded();
  return plot.evaluate((el,values)=>{
   const f=el._fullLayout,b=el.getBoundingClientRect();
   return {x0:b.x+f.xaxis._offset+f.xaxis.d2p(values[0]),x1:b.x+f.xaxis._offset+f.xaxis.d2p(values[1]),y:b.y+f.yaxis._offset+f.yaxis._length*.5};
  },[x0,x1]);
 };
 let p=await coords('.explorer .node-card .js-plotly-plot',19,26);
 await page.mouse.move(p.x0,p.y);await page.mouse.down();await page.mouse.move(p.x1,p.y,{steps:15});await page.mouse.up();
 await page.locator('.sel-handles').first().waitFor();
 await page.getByRole('button',{name:'X experiments',exact:true}).click();
 await page.getByRole('button',{name:'Use active window',exact:true}).click();
 const importedStart=await page.getByLabel('Experiment start bin',{exact:true}).inputValue();
 if (+importedStart<100) throw new Error('Selected window not mapped to chromosome bins');
 // Baseline manual flow begins unsplit and permits a non-argmax boundary.
 await page.getByRole('button',{name:'Manual segmentation',exact:true}).click();
 await page.getByRole('button',{name:'Place breakpoint',exact:true}).click();
 p=await coords('.manual-main .node-card .js-plotly-plot',20,20);
 await page.mouse.click(p.x0,p.y);
 await page.getByText(/manual split.*boundary 169/).waitFor();
 await page.locator('.split-status .badge.fail').waitFor({timeout:30000});
 await page.getByRole('button',{name:/undo \(1\)/}).click();
 await page.getByRole('button',{name:/redo \(1\)/}).click();
 await page.getByRole('button',{name:'Finish & assign CN (2 segments)',exact:true}).click();
 await page.locator('.card.cn').getByText(/reproduces exported CN row/).waitFor();
 await page.evaluate(()=>window.scrollTo(0,0));
 await page.screenshot({path:path.resolve(__dirname,'../validation/manual.png'),fullPage:true});
 // Compare estimators on unchanged manual boundaries and restore the option.
 const estimator=page.getByLabel('CN segment summary',{exact:true});
 for(const method of ['iqr_upper','iqr_two_sided']) {
  await estimator.selectOption(method);
  await page.locator('.card.cn').getByText(/OUTDATED/).waitFor();
  await page.getByRole('button',{name:'Finish & assign CN (2 segments)',exact:true}).click();
  await page.waitForFunction(()=>!document.querySelector('.cn.outdated') && document.querySelector('.estimator-comparison'));
  await page.locator('.card.cn').getByRole('columnheader',{name:'retained / total',exact:true}).waitFor();
 }
 const downloadPromise=page.waitForEvent('download');
 await page.getByRole('button',{name:'download JSON',exact:true}).click();
 const download=await downloadPromise, sessionPath=path.resolve(__dirname,'../validation/manual-iqr-session.json');
 await download.saveAs(sessionPath);
 const saved=JSON.parse(fs.readFileSync(sessionPath,'utf8')).state.manual['TTAGGCTAGGCCGGAA-1|chr9'];
 if(saved.cnEstimator!=='iqr_two_sided'||saved.cnResult.segment_estimator!=='iqr_two_sided'||saved.cnOutdated) throw new Error('Incorrect saved IQR state');
 const csvDownload=page.waitForEvent('download');
 await page.locator('.card.cn').getByRole('button',{name:'CSV',exact:true}).click();
 await (await csvDownload).saveAs(path.resolve(__dirname,'../validation/manual-iqr.csv'));
 await page.reload();
 await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.locator('input[type=file]').setInputFiles(sessionPath);
 await page.waitForFunction(()=>document.querySelector('select[aria-label="CN segment summary"]')?.value==='iqr_two_sided');
 await page.locator('.estimator-comparison').waitFor();
 await page.locator('.card.cn').scrollIntoViewIfNeeded();
 await page.screenshot({path:path.resolve(__dirname,'../validation/manual-iqr.png'),fullPage:false});
 // Drag directly on an experiment result to define the next source region.
 await page.locator('input[type=file]').setInputFiles(path.resolve(__dirname,'../validation/browser-session.json'));
 await page.getByRole('heading',{name:'X-value experiments',exact:true}).waitFor();
 const before=await page.getByLabel('Experiment start bin',{exact:true}).inputValue();
 p=await coords('.experiment-results .js-plotly-plot',70,80);
 await page.mouse.move(p.x0,p.y);await page.mouse.down();await page.mouse.move(p.x1,p.y,{steps:15});await page.mouse.up();
 await page.waitForFunction(old=>document.querySelector('input[aria-label="Experiment start bin"]').value!==old,before);
 if(errors.length)throw new Error(errors.join('\n'));
 fs.writeFileSync(path.resolve(__dirname,'../validation/workflow-results.json'),JSON.stringify({explorer_window_to_experiment:true,non_argmax_manual_split:169,undo_redo:true,full_genome_cn:true,iqr_both_methods:true,iqr_stale_state:true,iqr_session_restore:true,iqr_csv:true,result_region_drag:true,errors},null,2));
 console.log('PASS: explorer selection, experiment import, manual non-argmax split/tests/undo/redo/CN, result-region drawing');
 await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
