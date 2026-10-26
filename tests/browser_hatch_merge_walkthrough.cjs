// Merge-rules + Rule 0 browser walkthrough (2026-10-02, extended 2026-10-26): visualization/merge split, rule editor,
// live per-rule hatches, reorder, run, stale notice, session restore, legacy migration.
// Requires a running server (BASE_URL defaults to http://127.0.0.1:8766; the Rule 0 runs used a private server: python3 run.py --port 8767 with BASE_URL=http://127.0.0.1:8767) and Playwright + Chromium, e.g.:
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

 // ---- Rule 0 (pinned, built from the section-1 controls) ----
 const r0Parts=async()=>{const t=await panel.getByLabel('Rule 0 status').innerText();const m=t.match(/would merge now: (\d+)/),g=t.match(/blocked by gap: (\d+)/);return m?{m:+m[1],g:g?+g[1]:0}:null};
 // total Rule 0 would-merge segments = would merge now + blocked by gap
 const r0Flags=async()=>{const p=await r0Parts();return p?p.m+p.g:null};
 const tau=page.getByLabel('Rule 0 tau'),combine=page.getByLabel('Rule 0 combine');
 const srdCrit=page.getByLabel('Rule 0 criterion min |SRD/√φ| flank');
 const mergeSigns=page.getByLabel('Merge sign-disagreeing segments'),mergeFlanks=page.getByLabel('Merge flank-conflict segments');
 const waitR0Flags=()=>page.waitForFunction(()=>/would merge now: \d+/.test(document.querySelector('[aria-label="Rule 0 status"]')?.textContent||''),{timeout:30000});
 const restoreR0=async()=>{await tau.fill('0.2');await combine.selectOption('all')};
 const legendCounts=async()=>{const t=await page.locator('.rule0-legend').innerText();const w=t.match(/would merge now: (\d+)/),h=t.match(/held back: (\d+)/),g=t.match(/blocked by gap: (\d+)/);const wm=w?+w[1]:null,gap=g?+g[1]:0;return {wm,hb:h?+h[1]:null,gap,tot:wm===null?null:wm+gap}};
 const runEnabledWithRule0Only=!(await runBtn.isDisabled());

 // Rule 0 responds to its settings
 await combine.selectOption('any');await tau.fill('5');await waitR0Flags();await page.waitForTimeout(400);
 const r0High=await r0Flags();
 await tau.fill('0.0001');await srdCrit.uncheck();await page.waitForTimeout(500);
 const r0Low=await r0Flags();
 await restoreR0();await srdCrit.check();await page.waitForTimeout(400);

 // each merge switch on its own, then both: legend counts + the held-back reasons in the hover text
 const heldReasons=()=>plot.evaluate(el=>{const c={signs:0,flanks:0,both:0};for(const d of el.data){if(d.name!=='small segment')continue;const m=String(d.hovertemplate||'').match(/⚠ held back \(([^)]*)\)/);if(!m)continue;const sg=m[1].includes('signs disagree'),fl=m[1].includes('flanks disagree');if(sg&&fl)c.both++;else if(sg)c.signs++;else if(fl)c.flanks++}return c});
 const state=async(sg,fl)=>{await (sg?mergeSigns.check():mergeSigns.uncheck());await (fl?mergeFlanks.check():mergeFlanks.uncheck());await page.waitForTimeout(700);return {...(await legendCounts()),reasons:await heldReasons()}};
 // tau 5 + combine any only has sign holds; SRD/√φ ref 8 + combine all gives one sign-only and one sign+flank hold on chr8
 const srdRef=page.getByLabel('Rule 0 SRD/√φ ref');
 await srdRef.fill('8');await tau.fill('5');await combine.selectOption('all');
 const swOffOff=await state(false,false),swSignOnly=await state(true,false),swFlankOnly=await state(false,true),swBoth=await state(true,true);
 const {tot:wmOff,hb:hbOff}=swOffOff,{tot:wmOn,hb:hbOn}=swBoth;
 await state(false,false);
 // chr8 with SRD/√φ ref 8, τ 5, combine all (the 13-segment configuration; true defaults flag none on chr8), gap crossing off: every would-merge segment is gap-blocked, so a run merges 0 and nothing carries the x hatch
 const defaultsOff={...await r0Parts(),legend:await legendCounts(),notice:await page.getByText("All Rule 0 segments here are next to a genomic gap; tick 'allow merging across a genomic gap' in section 2 to merge them.").isVisible().catch(()=>false)};
 const allowGap=panel.getByLabel('Allow merging across a genomic gap');
 await allowGap.check();await page.waitForTimeout(800);
 const defaultsOn={...await r0Parts(),legend:await legendCounts(),notice:await page.getByText("All Rule 0 segments here are next to a genomic gap").isVisible().catch(()=>false)};
 await allowGap.uncheck();await page.waitForTimeout(500);
 await srdRef.fill('3.3');await restoreR0();await page.waitForTimeout(400);

 // combine alone: same tau and criteria, all vs any
 const combineTried=[];let combineAll=null,combineAny=null;
 for(const t of ['0.2','0.5','1','5']){await tau.fill(t);await combine.selectOption('all');await page.waitForTimeout(500);const a=await r0Flags();await combine.selectOption('any');await page.waitForTimeout(500);const b=await r0Flags();combineTried.push({tau:t,all:a,any:b});if(combineAll===null||b>a){combineAll=a;combineAny=b}if(b>a)break}
 await restoreR0();await page.waitForTimeout(300);

 // invalid Rule 0 blocks the run
 await tau.fill('');await page.waitForTimeout(400);
 const invalidR0Status=await panel.getByLabel('Rule 0 status').innerText();
 const runDisabledWhenR0Invalid=await runBtn.isDisabled();
 await tau.fill('0.2');await page.waitForTimeout(300);

 // non-pilot chromosome: live scores still give small segments
 const chromSel=tracks.locator('label').filter({hasText:/^chromosome/}).locator('select');
 await chromSel.selectOption('chr3');
 await page.waitForFunction(()=>String(document.querySelector('#cluster-tracks .js-plotly-plot')?.layout?.title?.text||'').includes('chr3'),{timeout:30000});
 await page.waitForTimeout(600);
 const chr3SmallSegments=await plot.evaluate(el=>el.data.filter(d=>d.name==='small segment').length);
 await chromSel.selectOption('chr8');
 await page.waitForFunction(()=>String(document.querySelector('#cluster-tracks .js-plotly-plot')?.layout?.title?.text||'').includes('chr8'),{timeout:30000});
 await page.waitForTimeout(600);

 // candidate audit shows derived p-value reasons
 const auditBox=tracks.getByLabel('candidate audit (rejected / blocked)');
 await auditBox.check();await page.waitForTimeout(800);
 const auditDerivedReason=await plot.evaluate(el=>el.data.some(d=>[].concat(d.hovertemplate||[]).some(h=>String(h).includes('permutation p ='))));
 await auditBox.uncheck();await page.waitForTimeout(300);

 // a Rule 0-only run names Rule 0 in the history
 // chr8's would-merge segments all sit next to genomic gaps, so a run without gap crossing makes 0 merges
 await combine.selectOption('any');await tau.fill('5');await panel.getByLabel('Allow merging across a genomic gap').check();await page.waitForTimeout(400);
 await runBtn.click();
 await panel.locator('table.tbl, p[role="status"]').first().waitFor({timeout:60000});
 await page.waitForTimeout(1500);
 const rule0HistoryRows=await panel.locator('table.tbl tbody tr').count();
 const rule0FirstRuleCell=rule0HistoryRows?await panel.locator('table.tbl tbody tr').first().locator('td').nth(6).innerText():null;
 await restoreR0();await page.waitForTimeout(300);

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
 const rulesAfterRestore=await panel.locator('.rule-card:not(.rule0-card)').count();

 // legacy session: only hatchLayers present -> three migrated preset rules
 const session=JSON.parse(fs.readFileSync(sessionPath,'utf8'));
 const owner=findOwner(session,'mergeRules');delete owner.mergeRules;delete owner.allowGapCrossing;
 owner.seg={...owner.seg,metric:'log2fc',flank:'min_srd',labels:{...owner.seg.labels,proposal:true,consistency:false}};delete owner.seg.labels.signs;delete owner.seg.colorMetric;delete owner.seg.criteria;delete owner.seg.signs;
 owner.hatchLayers={transition:{on:true,metric:'srd_phi',estimator:'mean'},proposal:{on:false,metric:'srd_phi',estimator:'mean',threshold:2.5},ambiguous:{on:true,estimator:'median'}};
 const legacyPath=path.join(dir,'legacy.json');fs.writeFileSync(legacyPath,JSON.stringify(session));
 await fileInput.setInputFiles(legacyPath);await page.waitForTimeout(600);
 const nCards=await panel.locator('.rule-card:not(.rule0-card)').count();const legacyNames=[],legacyEnabled=[];
 for(let i=1;i<=nCards;i++){legacyNames.push(await panel.getByLabel(`Rule ${i} name`).inputValue());legacyEnabled.push(await panel.getByLabel(`Rule ${i} enabled`).isChecked())}
 const colourByAfterLegacy=await tracks.locator('label').filter({hasText:/^colour by/}).locator('select').inputValue();
 const signsLayerAfterLegacy=await tracks.getByLabel('signs-disagree hatch').isChecked();
 fs.rmSync(dir,{recursive:true,force:true});

 const result={generated_at:new Date().toISOString(),visHeading,shapesBaseline,shapesWithRule,runEnabledWithRule0Only,r0High,r0Low,wmOff,hbOff,wmOn,hbOn,defaultsOff,defaultsOn,swOffOff,swSignOnly,swFlankOnly,swBoth,combineAll,combineAny,combineTried,invalidR0Status,runDisabledWhenR0Invalid,chr3SmallSegments,auditDerivedReason,rule0HistoryRows,rule0FirstRuleCell,colourByAfterLegacy,signsLayerAfterLegacy,flagsHigh,flagsLow,emptyThresholdProblem,runDisabledWhenIncomplete,legendHasBothRules,rule1NameAfterMove,summaryText,headers,historyRows,firstRuleCell,staleText,enabledBeforeRestore,enabledAfterRestore,rulesAfterRestore,legacyNames,legacyEnabled,errors};
 fs.writeFileSync(path.join(out,'merge_rules_check_results.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result,null,2));
 assert.ok(visHeading>0,'visualization section heading missing');
 assert.equal(runEnabledWithRule0Only,true,'Rule 0 alone should allow a run');
 assert.ok(r0High>r0Low&&r0Low===0,'Rule 0 should flag more at tau 5 than at tau 0.0001 without the SRD criterion');
 assert.ok(swOffOff.reasons.both>0&&swOffOff.reasons.signs>0,'test data should hold segments for signs only and for signs+flanks');
 assert.ok(swOffOff.tot>0&&swOffOff.gap>0,'default (gap crossing off) legend should show blocked-by-gap');
 assert.ok(defaultsOff.m===0&&defaultsOff.g>0,'chr8 at defaults with gap crossing off: would merge now 0, blocked by gap > 0');
 assert.equal(defaultsOff.legend.wm,0,'legend agrees: would merge now 0');assert.equal(defaultsOff.legend.gap,defaultsOff.g,'legend and Rule 0 card agree on blocked by gap');
 assert.equal(defaultsOff.notice,true,'the all-gap-blocked notice should be visible');
 assert.ok(defaultsOn.m>0&&defaultsOn.g===0,'allowing gap crossing: would merge now > 0, blocked by gap 0');
 assert.equal(defaultsOn.m,defaultsOff.g,'allowing gap crossing releases exactly the gap-blocked segments');
 assert.equal(defaultsOn.notice,false,'the notice disappears once gap crossing is allowed');
 assert.ok(swSignOnly.hb<=swOffOff.hb&&swFlankOnly.hb<=swOffOff.hb,'a single switch must not increase held-back');
 assert.equal(swSignOnly.reasons.signs+swSignOnly.reasons.flanks+swSignOnly.reasons.both,swSignOnly.hb,'hover held-back reasons should match the legend');
 assert.equal(swSignOnly.hb,swOffOff.reasons.flanks+swOffOff.reasons.both,'sign switch on: only flank-conflict holds remain');
 assert.equal(swFlankOnly.hb,swOffOff.reasons.signs+swOffOff.reasons.both,'flank switch on: only sign holds remain');
 assert.ok(combineAny>=combineAll,'combine any should flag at least as many segments as all');
 assert.equal(hbOn,0,'with both merge switches on nothing is held back');
 assert.equal(wmOn,wmOff+hbOff,'switching both on should move every held-back segment to would-merge');
 assert.ok(invalidR0Status.includes('τ'),'a cleared tau should be reported in the Rule 0 status');
 assert.equal(runDisabledWhenR0Invalid,true,'invalid Rule 0 should block the run');
 assert.ok(chr3SmallSegments>0,'chr3 should show small segments from live scores');
 assert.equal(auditDerivedReason,true,'candidate audit hovers should include derived permutation p reasons');
 assert.ok(rule0HistoryRows>0,'a Rule 0-only run with gap crossing should merge something');
 assert.ok(rule0FirstRuleCell.includes('Rule 0'),'Rule 0 merges should be named in the history');
 assert.equal(colourByAfterLegacy,'delta_mean','legacy seg metric should migrate to delta_mean');
 assert.equal(signsLayerAfterLegacy,false,'legacy consistency=false should leave the signs layer off');
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
