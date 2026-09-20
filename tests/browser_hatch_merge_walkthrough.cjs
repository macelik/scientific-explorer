// Task 10 (2026-09-20) browser walkthrough for the exploratory hatch layers and the
// adjacent-segment merge panel (HatchMergeView.tsx / the "new exploratory hatch layers"
// controls in IntegrationSegmentationView.tsx). Requires a running server
// (python3 run.py --port 8766) and Playwright + Chromium, e.g.:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright \
//   CHROMIUM_PATH=/path/to/chrome \
//   node tests/browser_hatch_merge_walkthrough.cjs
// Writes validation/hatch-merge/hatch_merge_check_results.json and
// validation/hatch-merge/merge-step10-hatch.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path'),os=require('os'),assert=require('assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:1000}}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console', msg => { if (msg.type() === 'error') errors.push('console: '+msg.text()) });
 await page.goto('http://127.0.0.1:8766');await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});

 // ---- Chromosome tracks: baseline shape count, then each new hatch layer individually ----
 const tracks=page.locator('#cluster-tracks');await tracks.scrollIntoViewIfNeeded();
 const plot=tracks.locator('.js-plotly-plot');await plot.waitFor();
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?._fullLayout);
 const legend = () => page.locator('.seg-legend').last().innerText();
 const shapeCount = () => plot.evaluate(el=>el.layout.shapes.length);

 const shapesBaseline = await shapeCount();

 await page.getByRole('checkbox',{name:'transition zone (/)',exact:true}).check();
 await page.waitForTimeout(600);
 const shapesAfterTransition = await shapeCount();
 const legendAfterTransition = await legend();

 await page.getByRole('checkbox',{name:'merge proposal (×)',exact:true}).check();
 await page.waitForTimeout(600);
 const shapesAfterProposal = await shapeCount();
 const legendAfterProposal = await legend();

 await page.getByRole('checkbox',{name:'ambiguous flank preference (×)',exact:true}).check();
 await page.waitForTimeout(600);
 const shapesAfterAmbiguous = await shapeCount();
 const legendAfterAmbiguous = await legend();

 const hatchToggles = {
   shapesBaseline, shapesAfterTransition, shapesAfterProposal, shapesAfterAmbiguous,
   legendMentionsTransition: legendAfterTransition.includes('opposite signs'),
   legendMentionsProposal: legendAfterProposal.includes("at least one flank's"),
   legendMentionsAmbiguous: legendAfterAmbiguous.includes('differs from the weaker'),
 };

 // ---- Exploratory adjacent-segment merge panel ----
 const mergePanel = page.locator('.hatch-merge');
 await mergePanel.scrollIntoViewIfNeeded();
 await page.getByRole('heading',{name:'Exploratory adjacent-segment merge'}).waitFor();
 const thresholdInput = mergePanel.locator('input[aria-label="Merge threshold"]');
 await thresholdInput.fill('50');
 await mergePanel.getByRole('checkbox',{name:'allow merging across a genomic gap (off by default)'}).check();
 await page.getByRole('button',{name:'Run merge'}).click();
 await page.waitForFunction(() => document.querySelectorAll('.hatch-merge .js-plotly-plot').length >= 2, {timeout: 30000});
 await page.waitForTimeout(500);

 const plots = mergePanel.locator('.js-plotly-plot');
 const nPlots = await plots.count();
 const summaryText = await mergePanel.locator('p[role="status"]').innerText();

 // Move the step slider to the final step so the merged-interval hatch overlay (the
 // fixed yref) is present, and inspect the second plot's rendered shapes directly.
 const slider = mergePanel.locator('input[aria-label="Merge step"]');
 const maxStep = await slider.getAttribute('max');
 await slider.fill(String(maxStep));
 await slider.dispatchEvent('input');
 await page.waitForTimeout(800);

 const secondPlotShapes = await plots.nth(1).evaluate(el => el.layout.shapes.map(s => ({type:s.type, yref:s.yref})));
 const hatchLineShapes = secondPlotShapes.filter(s => s.type === 'line');
 const badYref = hatchLineShapes.some(s => typeof s.yref === 'string' && s.yref.trim().startsWith('domain'));
 const stepLabel = await mergePanel.locator('.muted.small').first().innerText();

 const jsonBtn = await mergePanel.getByRole('button', {name: 'Download JSON'}).count();
 const csvBtn = await mergePanel.getByRole('button', {name: 'Download merge-history CSV'}).count();

 const out=path.resolve(__dirname, '../validation/hatch-merge');
 fs.mkdirSync(out,{recursive:true});
 await mergePanel.screenshot({path: path.join(out, 'merge-step10-hatch.png')});

 const historyRows = await mergePanel.locator('table.tbl tbody tr').count();

 // Reset / Undo
 await mergePanel.getByRole('button',{name:'Reset to original'}).click();
 await page.waitForTimeout(300);
 const afterReset = await mergePanel.locator('.muted.small').first().innerText();
 await mergePanel.getByRole('button',{name:'Jump to final'}).click();
 await page.waitForTimeout(300);
 await mergePanel.getByRole('button',{name:'Undo one step'}).click();
 await page.waitForTimeout(300);
 const afterUndo = await mergePanel.locator('.muted.small').first().innerText();

 // Stale-result notice: change threshold after a run without rerunning
 await thresholdInput.fill('10');
 await page.waitForTimeout(300);
 const staleText = await mergePanel.locator('p[role="status"]').innerText();
 const staleClass = await mergePanel.locator('p[role="status"]').getAttribute('class');

 // Chromosome-change responsiveness: switch the shared chromosome selector and confirm
 // the tracks plot updates without hanging.
 const chromSelect = tracks.locator('select').first();
 const chromBefore = await chromSelect.inputValue();
 await chromSelect.selectOption('chr9');
 await page.waitForFunction(() => {
   const el = document.querySelector('#cluster-tracks .js-plotly-plot');
   return el && el._fullLayout && el._fullLayout.title && el._fullLayout.title.text && el._fullLayout.title.text.includes('chr9');
 }, {timeout: 20000});
 const chromAfter = await chromSelect.inputValue();
 await chromSelect.selectOption('chr8');

 // Session restore: hatchLayers.transition.on should round-trip via download JSON -> upload JSON.
 // Field names are the literal checkbox-checked state at each point (not double-negated).
 const transitionCb = page.getByRole('checkbox',{name:'transition zone (/)',exact:true});
 await transitionCb.check();
 const checkedBeforeExport = await transitionCb.isChecked();
 const [download] = await Promise.all([
   page.waitForEvent('download'),
   page.getByRole('button',{name:'download JSON',exact:true}).click(),
 ]);
 const sessionPath = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'hatch-merge-session-')), 'session-export.json');
 await download.saveAs(sessionPath);
 await transitionCb.uncheck();
 const checkedBeforeRestore = await transitionCb.isChecked(); // expected false
 const fileInput = page.locator('input[type="file"][accept="application/json"]');
 await fileInput.setInputFiles(sessionPath);
 await page.waitForTimeout(500);
 const checkedAfterRestore = await transitionCb.isChecked(); // expected true if restore worked
 fs.rmSync(path.dirname(sessionPath), {recursive: true, force: true});

 const result = {
   generated_at: new Date().toISOString(),
   hatchToggles,
   mergeRun: { summaryText, nPlots, secondPlotShapeCount: secondPlotShapes.length,
     hatchLineShapeCount: hatchLineShapes.length, hatchLineYrefSample: hatchLineShapes.slice(0,3),
     badYrefFound: badYref, stepLabel, jsonBtn, csvBtn, historyRows },
   resetUndo: { afterReset, afterUndo },
   staleNotice: { staleText, staleClass },
   chromosomeChange: { chromBefore, chromAfter },
   sessionRestore: { checkedBeforeExport, checkedBeforeRestore, checkedAfterRestore },
   errors,
 };
 fs.writeFileSync(path.join(out,'hatch_merge_check_results.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result,null,2));
 assert.equal(badYref, false, 'Found an invalid " domain" yref in the merged-interval hatch shapes');
 assert.equal(checkedBeforeRestore, false, 'checkbox should read unchecked right before the restore');
 assert.equal(checkedAfterRestore, true, 'session restore did not re-check the hatch layer');
 assert.deepEqual(errors, [], 'console/page errors present');
 console.log('PASS hatch layers + merge panel + yref fix + session restore verified');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
