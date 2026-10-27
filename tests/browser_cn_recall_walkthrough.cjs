// Genome merge + CN re-call walkthrough (2026-10-27): button, status, Before/After/Changed on both karyograms,
// impact summary, export, stale notice.
// Requires a private server (python3 run.py --port 8799; BASE_URL=http://127.0.0.1:8799) and Playwright + Chromium:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome BASE_URL=http://127.0.0.1:8799 \
//   node tests/browser_cn_recall_walkthrough.cjs
// Writes validation/cn-recall/cn_recall_check_results.json and validation/cn-recall/cn-recall.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const BASE_URL=process.env.BASE_URL||'http://127.0.0.1:8766';
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
const PALETTE=[[123,31,162],[206,147,216],[76,175,80],[255,138,101],[171,0,0]];
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:1100},acceptDownloads:true}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console',msg=>{if(msg.type()==='error')errors.push('console: '+msg.text())});
 await page.goto(BASE_URL);await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});
 const clusterPanel=page.locator('.panel',{has:page.getByRole('heading',{name:/Cluster-level CN/})});
 const clusterCanvas=clusterPanel.locator('canvas').first();
 const cellPanel=page.locator('.panel',{has:page.getByRole('heading',{name:/Per-cell CN under the cluster segmentation/})});
 // pixels exactly equal to a palette colour inside the plot area (labels/bands excluded: x >= 164, y >= 20)
 const palettePixels=()=>clusterCanvas.evaluate((cv,P)=>{const d=cv.getContext('2d').getImageData(164,20,cv.width-164,cv.height-20).data;let n=0;for(let i=0;i<d.length;i+=4)if(P.some(p=>p[0]===d[i]&&p[1]===d[i+1]&&p[2]===d[i+2]))n++;return n},PALETTE);
 const digest=()=>clusterCanvas.evaluate(cv=>{const d=cv.getContext('2d').getImageData(0,0,cv.width,cv.height).data;let h=0;for(let i=0;i<d.length;i+=7)h=(h*31+d[i])>>>0;return h});
 const modeBtn=n=>page.getByRole('group',{name:'Karyogram mode'}).getByRole('button',{name:n,exact:true});
 const afterDisabledInitially=await modeBtn('After').isDisabled(),changedDisabledInitially=await modeBtn('Changed').isDisabled();
 const beforeDigest=await digest(),beforePalette=await palettePixels();

 await page.getByLabel('Rule 0 tau').fill('0.5');await page.getByLabel('Rule 0 SRD/√φ ref').fill('8');await page.getByLabel('Rule 0 combine').selectOption('any');
 const panel=page.locator('.hatch-merge');await panel.scrollIntoViewIfNeeded();
 const t0=Date.now();
 await panel.getByRole('button',{name:'Merge genome and re-call CN'}).click();
 const status=page.getByLabel('Genome re-call status');
 await page.waitForFunction(()=>/cluster CN ✓|error:/.test(document.querySelector('[aria-label="Genome re-call status"]')?.textContent||''),undefined,{timeout:300000});
 const midText=await status.innerText();
 const afterEnabledMid=!(await modeBtn('After').isDisabled());
 await modeBtn('After').click();
 const perCellPendingShown=(await cellPanel.getByText(/re-calling per-cell CN/).count())>0||/per-cell ✓/.test(await status.innerText());
 await page.waitForFunction(()=>/per-cell ✓ · done|error:/.test(document.querySelector('[aria-label="Genome re-call status"]')?.textContent||''),undefined,{timeout:600000});
 const doneText=await status.innerText(),elapsedS=Math.round((Date.now()-t0)/1000);

 const impact=page.getByLabel('CN re-call impact');
 const rows=await impact.locator('tbody tr').evaluateAll(trs=>trs.map(tr=>[...tr.children].map(td=>td.textContent)));
 const merges=rows.reduce((s,r)=>s+parseInt(r[2],10),0),binsChanged=rows.reduce((s,r)=>s+parseInt(r[4].replace(/,/g,''),10),0);
 const cellLine=await impact.getByText(/cells with any change:/).innerText();
 await cellPanel.locator('canvas').first().waitFor({timeout:60000}).catch(()=>{});// percell.bin is fetched after the status flips to done
 const perCellCanvasAfter=await cellPanel.locator('canvas').count();
 await modeBtn('Changed').click();await page.waitForTimeout(400);
 const changedPalette=await palettePixels();
 await clusterPanel.screenshot({path:path.join('validation','cn-recall','cn-recall.png')}).catch(async()=>{fs.mkdirSync(path.join('validation','cn-recall'),{recursive:true});await clusterPanel.screenshot({path:path.join('validation','cn-recall','cn-recall.png')})});
 const [dl]=await Promise.all([page.waitForEvent('download'),impact.getByRole('button',{name:'merged_breakpoints.tsv'}).click()]);
 const exported=fs.readFileSync(await dl.path(),'utf8').split('\n');
 await modeBtn('Before').click();await page.waitForTimeout(400);
 const backDigest=await digest();
 await page.getByLabel('Rule 0 tau').fill('0.4');await page.waitForTimeout(300);
 const staleAfterEdit=(await status.innerText()).includes('Rules changed since this re-call; showing the previous one.');
 // ---- merge rounds: adopt -> tracks lose breakpoints -> noise table gains Round 1 -> same rules give 0 merges -> undo ----
 const bpCount=async()=>{await page.locator('#cluster-tracks .deferred-plot').first().scrollIntoViewIfNeeded();await page.locator('#cluster-tracks .js-plotly-plot').waitFor({timeout:60000});await page.waitForTimeout(500);return bpCountNow()};
 const bpCountNow=()=>page.evaluate(()=>{const el=document.querySelector('#cluster-tracks .js-plotly-plot');return ((el&&el.layout&&el.layout.shapes)||[]).filter(sh=>sh.type==='line'&&String(sh.yref||'').endsWith(' domain')).length});
 const baseStrip=page.getByLabel('Segmentation base');
 const noise=page.getByLabel('Call noise by round');
 // value of a measure row in the first source block's cluster level, column named colName (counted from the right: first rows carry extra rowSpan cells)
 const noiseValue=(colName)=>noise.locator('table').first().evaluate((t,colName)=>{
   const heads=[...t.querySelectorAll('thead th')].map(h=>h.textContent.trim());const hi=heads.indexOf(colName);if(hi<0)return null;
   const fromRight=heads.length-hi;
   for(const tr of t.querySelectorAll('tbody tr')){const cells=[...tr.children].map(c=>c.textContent.trim());
     if(cells.includes('zig-zags'))return parseInt(cells[cells.length-fromRight].replace(/,/g,''),10)}
   return null},colName);
 const bpBefore=await bpCount();
 await noise.scrollIntoViewIfNeeded();
 await noise.getByRole('group',{name:'Noise scope'}).getByRole('button',{name:'genome-wide',exact:true}).click();
 await page.waitForTimeout(300);
 const zzOriginal=await noiseValue('Original');
 await page.getByLabel('Rule 0 tau').fill('0.5');await page.getByLabel('Rule 0 SRD/√φ ref').fill('8');await page.getByLabel('Rule 0 combine').selectOption('any');
 await impact.getByRole('button',{name:'Adopt as round 1'}).click();
 await page.waitForFunction(()=>/Round 1 \(/.test(document.querySelector('[aria-label="Segmentation base"]')?.textContent||''),undefined,{timeout:120000});
 await page.waitForFunction(()=>[...document.querySelectorAll('[aria-label="Call noise by round"] thead th')].some(h=>h.textContent.trim()==='Round 1'),undefined,{timeout:120000});
 await page.waitForTimeout(1500);
 const bpAfter=await bpCount();
 const zzRound1=await noiseValue('Round 1');
 const defs=page.getByLabel('Noise definitions');
 const definitionsVisible=await defs.isVisible()&&(await defs.innerText()).includes('Zig-zag: an internal segment');
 // re-call on the base with the same rules: every round merges on the original pseudobulk, so nothing new merges
 await panel.scrollIntoViewIfNeeded();
 await panel.getByRole('button',{name:'Merge genome and re-call CN'}).click();
 await page.waitForFunction(()=>/per-cell ✓ · done|error:/.test(document.querySelector('[aria-label="Genome re-call status"]')?.textContent||''),undefined,{timeout:600000});
 const status2=await status.innerText();
 const rows2=await page.getByLabel('CN re-call impact').locator('tbody tr').evaluateAll(trs=>trs.map(tr=>[...tr.children].map(td=>td.textContent)));
 const round2Merges=rows2.reduce((s,r)=>s+parseInt(r[2],10),0);
 await page.getByRole('button',{name:'Undo last round'}).click();
 await page.waitForFunction(()=>(document.querySelector('[aria-label="Segmentation base"]')?.textContent||'').trim()==='Segmentation base: Original',undefined,{timeout:120000});
 await page.waitForTimeout(1500);
 const bpUndo=await bpCount();
 const rounds={bpBefore,bpAfter,bpUndo,zzOriginal,zzRound1,round2Merges,definitionsVisible};
 await page.getByLabel('Rule 0 tau').fill('0.2');await page.getByLabel('Rule 0 SRD/√φ ref').fill('3.3');await page.getByLabel('Rule 0 combine').selectOption('all');

 const result={date:new Date().toISOString(),settings:{tau:0.5,srdRef:8,combine:'any'},elapsedS,midText,doneText,merges,binsChanged,cellLine,
   afterDisabledInitially,changedDisabledInitially,afterEnabledMid,perCellPendingShown,perCellCanvasAfter,beforePalette,changedPalette,
   backToBeforeIdentical:backDigest===beforeDigest,exportHeader:exported.slice(0,4),staleAfterEdit,rounds,errors};
 fs.mkdirSync(path.join('validation','cn-recall'),{recursive:true});
 fs.writeFileSync(path.join('validation','cn-recall','cn_recall_check_results.json'),JSON.stringify(result,null,2));
 assert.equal(afterDisabledInitially,true);assert.equal(changedDisabledInitially,true);
 assert.ok(!doneText.includes('error:'),doneText);
 assert.equal(afterEnabledMid,true,'After enables once cluster CN is ready');
 assert.equal(perCellPendingShown,true);
 assert.ok(merges>0,'a loose Rule 0 should merge something genome-wide');
 assert.ok(binsChanged>0,'re-calling should change at least one cluster bin');
 assert.ok(perCellCanvasAfter>=1,'per-cell karyogram should draw in After mode');
 assert.ok(changedPalette>0&&changedPalette<beforePalette,'Changed draws only changed bins in colour');
 assert.equal(result.backToBeforeIdentical,true,'Before must equal the original rendering');
 assert.ok(exported[0].startsWith('# rules: ')&&exported[3].includes('absolute_bin'),'export header + columns');
 assert.equal(staleAfterEdit,true,'editing a rule after a re-call shows the stale notice');
 assert.ok(!status2.includes('error:'),status2);
 assert.ok(bpBefore>0&&bpAfter<bpBefore,'adopting removes breakpoints from the tracks: '+bpBefore+' -> '+bpAfter);
 assert.ok(zzOriginal!=null&&zzRound1!=null&&zzRound1<zzOriginal,'Round 1 has fewer genome-wide zig-zags: '+zzOriginal+' -> '+zzRound1);
 assert.equal(definitionsVisible,true,'noise definitions visible');
 assert.equal(round2Merges,0,'same rules on the base merge nothing');
 assert.equal(bpUndo,bpBefore,'undo restores the original breakpoints');
 assert.deepEqual(errors,[],'console/page errors present');
 console.log('PASS cn re-call: button, status, Before/After/Changed, summary, export, stale notice, merge rounds');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
