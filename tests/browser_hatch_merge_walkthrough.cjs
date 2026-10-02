// Merge-rules browser walkthrough (2026-10-02): visualization/merge split, rule editor,
// live per-rule hatches, reorder, run, stale notice, session restore, legacy migration.
// Requires a running server (python3 run.py --port 8766; override the page URL with BASE_URL) and Playwright + Chromium, e.g.:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome \
//   node tests/browser_hatch_merge_walkthrough.cjs
// Writes validation/hatch-merge/merge_rules_check_results.json and validation/hatch-merge/merge-rules.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const BASE_URL=process.env.BASE_URL||'http://127.0.0.1:8766';
const fs=require('fs'),path=require('path'),os=require('os'),assert=require('assert/strict');
function findOwner(obj,key){if(!obj||typeof obj!=='object')return null;if(Object.prototype.hasOwnProperty.call(obj,key))return obj;for(const v of Object.values(obj)){const f=findOwner(v,key);if(f)return f}return null}
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:1100}}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console',msg=>{if(msg.type()==='error')errors.push('console: '+msg.text())});
 await page.goto(BASE_URL);await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});

 const tracks=page.locator('#cluster-tracks');await tracks.scrollIntoViewIfNeeded();
 const plot=tracks.locator('.js-plotly-plot');await plot.waitFor();
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?._fullLayout);
 const shapeCount=()=>plot.evaluate(el=>el.layout.shapes.length);
 const visHeading=await tracks.getByText('1. Visualization').count();
 const shapesBaseline=await shapeCount();

 const panel=page.locator('.hatch-merge');await panel.scrollIntoViewIfNeeded();
 await page.getByRole('heading',{name:/2\. Merge rules/}).waitFor();
 const runBtn=panel.getByRole('button',{name:'Run merge'});
 const runDisabledWithNoRules=await runBtn.isDisabled();

 // rule 1: weak-flank preset; drive the threshold high, then near zero
 await panel.getByLabel('Add rule').selectOption('weak_flank');
 const status1=panel.getByLabel('Rule 1 status');
 const flags=async(label)=>{const m=(await panel.getByLabel(label).innerText()).match(/flags (\d+)/);return m?+m[1]:null};
 const thr1=panel.getByLabel('Rule 1 condition 1 threshold');
 await thr1.fill('1000');
 await page.waitForFunction(()=>/flags [1-9]\d*/.test(document.querySelector('[aria-label="Rule 1 status"]')?.textContent||''),{timeout:30000});
 const flagsHigh=await flags('Rule 1 status');
 await page.waitForTimeout(600);
 const shapesWithRule=await shapeCount();
 await thr1.fill('0.0001');await page.waitForTimeout(400);
 const flagsLow=await flags('Rule 1 status');
 await thr1.fill('');await page.waitForTimeout(300);
 const emptyThresholdProblem=(await status1.innerText()).includes('threshold');
 const runDisabledWhenIncomplete=await runBtn.isDisabled();
 await thr1.fill('3.3');await page.waitForTimeout(400);

 // rule 2: sign-mismatch preset; legend lists both; move it to the top
 await panel.getByLabel('Add rule').selectOption('sign_mismatch');await page.waitForTimeout(400);
 const legendText=await tracks.locator('.seg-legend').last().innerText();
 const legendHasBothRules=legendText.includes('weak SRD/√φ flank')&&legendText.includes('flank signs disagree');
 await panel.getByLabel('Move rule 2 up').click();
 const rule1NameAfterMove=await panel.getByLabel('Rule 1 name').inputValue();

 // run
 await panel.getByLabel('Allow merging across a genomic gap').check();
 await runBtn.click();
 // result plots are lazy (DeferredPlot): bring the result table, just below them, into view
 await panel.locator('table.tbl').waitFor();await panel.locator('table.tbl').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelectorAll('.hatch-merge .js-plotly-plot').length>=2,{timeout:60000});
 await page.waitForTimeout(500);
 const summaryText=await panel.locator('p[role="status"]').first().innerText();
 const headers=await panel.locator('table.tbl thead th').allInnerTexts();
 const historyRows=await panel.locator('table.tbl tbody tr').count();
 const firstRuleCell=historyRows?await panel.locator('table.tbl tbody tr').first().locator('td').nth(6).innerText():null;
 const out=path.resolve(__dirname,'../validation/hatch-merge');fs.mkdirSync(out,{recursive:true});
 await panel.screenshot({path:path.join(out,'merge-rules.png')});

 // stale notice after editing a rule (the weak-flank rule is now rule 2)
 await panel.getByLabel('Rule 2 condition 1 threshold').fill('2.0');await page.waitForTimeout(300);
 const staleText=await panel.locator('p[role="status"]').first().innerText();

 // session export -> change -> restore
 const enabled1=panel.getByLabel('Rule 1 enabled');
 const [download]=await Promise.all([page.waitForEvent('download'),page.getByRole('button',{name:'download JSON',exact:true}).click()]);
 const dir=fs.mkdtempSync(path.join(os.tmpdir(),'merge-rules-session-'));const sessionPath=path.join(dir,'session.json');
 await download.saveAs(sessionPath);
 await enabled1.uncheck();const enabledBeforeRestore=await enabled1.isChecked();
 const fileInput=page.locator('input[type="file"][accept="application/json"]');
 await fileInput.setInputFiles(sessionPath);await page.waitForTimeout(600);
 const enabledAfterRestore=await panel.getByLabel('Rule 1 enabled').isChecked();
 const rulesAfterRestore=await panel.locator('.rule-card').count();

 // legacy session: only hatchLayers present -> three migrated preset rules
 const session=JSON.parse(fs.readFileSync(sessionPath,'utf8'));
 const owner=findOwner(session,'mergeRules');delete owner.mergeRules;delete owner.allowGapCrossing;
 owner.hatchLayers={transition:{on:true,metric:'srd_phi',estimator:'mean'},proposal:{on:false,metric:'srd_phi',estimator:'mean',threshold:2.5},ambiguous:{on:true,estimator:'median'}};
 const legacyPath=path.join(dir,'legacy.json');fs.writeFileSync(legacyPath,JSON.stringify(session));
 await fileInput.setInputFiles(legacyPath);await page.waitForTimeout(600);
 const nCards=await panel.locator('.rule-card').count();const legacyNames=[],legacyEnabled=[];
 for(let i=1;i<=nCards;i++){legacyNames.push(await panel.getByLabel(`Rule ${i} name`).inputValue());legacyEnabled.push(await panel.getByLabel(`Rule ${i} enabled`).isChecked())}
 fs.rmSync(dir,{recursive:true,force:true});

 const result={generated_at:new Date().toISOString(),visHeading,shapesBaseline,shapesWithRule,runDisabledWithNoRules,flagsHigh,flagsLow,emptyThresholdProblem,runDisabledWhenIncomplete,legendHasBothRules,rule1NameAfterMove,summaryText,headers,historyRows,firstRuleCell,staleText,enabledBeforeRestore,enabledAfterRestore,rulesAfterRestore,legacyNames,legacyEnabled,errors};
 fs.writeFileSync(path.join(out,'merge_rules_check_results.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result,null,2));
 assert.ok(visHeading>0,'visualization section heading missing');
 assert.equal(runDisabledWithNoRules,true,'Run merge should be disabled with no rules');
 assert.ok(flagsHigh>flagsLow,'a high threshold should flag more segments than a near-zero one');
 assert.ok(shapesWithRule>shapesBaseline,'an enabled rule should add hatch shapes to the tracks');
 assert.equal(emptyThresholdProblem,true,'a cleared threshold should show a problem on the rule card');
 assert.equal(runDisabledWhenIncomplete,true,'an incomplete enabled rule should block the run');
 assert.equal(legendHasBothRules,true,'legend should list both rules');
 assert.equal(rule1NameAfterMove,'flank signs disagree','moving rule 2 up should put it first');
 assert.ok(headers.includes('rule'),'history table needs a rule column');
 if(historyRows>0)assert.ok(firstRuleCell&&firstRuleCell.trim().length>0,'history rows should name the rule');
 assert.ok(staleText.includes('Settings changed'),'editing a rule after a run should mark the result stale');
 assert.equal(enabledBeforeRestore,false);
 assert.equal(enabledAfterRestore,true,'session restore should bring back the saved rule state');
 assert.equal(rulesAfterRestore,2);
 assert.deepEqual(legacyNames,['transition zone (from saved session)','merge proposal (from saved session)','ambiguous flank preference (from saved session)'],'legacy hatchLayers should migrate to three preset rules');
 assert.deepEqual(legacyEnabled,[true,false,true]);
 assert.deepEqual(errors,[],'console/page errors present');
 console.log('PASS merge rules: split sections, editor, live hatches, reorder, run, stale, session restore, legacy migration');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
