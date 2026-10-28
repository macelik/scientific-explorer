// Allelic segmentation walkthrough (2026-10-28): tab, tracks before a run, |dBAF| run, tables, signal switch, depth original.
// Requires a private server (python3 run.py --port 8799; BASE_URL=http://127.0.0.1:8799) and Playwright + Chromium:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome BASE_URL=http://127.0.0.1:8799 \
//   node tests/browser_allelic_walkthrough.cjs
// Writes validation/allelic/allelic_check_results.json and validation/allelic/allelic.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const BASE_URL=process.env.BASE_URL||'http://127.0.0.1:8766';
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
  const page=await browser.newPage({viewport:{width:1400,height:1300},acceptDownloads:true}),errors=[];
  page.on('pageerror',e=>errors.push(String(e)));
  page.on('console',m=>{if(m.type()==='error')errors.push('console: '+m.text())});
  await page.goto(BASE_URL);await page.getByText('300 matching cells',{exact:true}).waitFor();
  await page.getByRole('button',{name:'Integration · allelic segmentation',exact:true}).click();
  await page.getByText('Run segmentation to see allelic segments.').waitFor({timeout:180000});
  await page.waitForFunction(()=>{const e=document.querySelector('.allelic-view .js-plotly-plot');return !!(e&&e.layout&&Object.keys(e.layout).filter(k=>/^yaxis\d*$/.test(k)).length>=10)},null,{timeout:120000});
  const panels=await page.evaluate(()=>Object.keys(document.querySelector('.allelic-view .js-plotly-plot').layout).filter(k=>/^yaxis\d*$/.test(k)).length);
  const t0=Date.now();
  await page.getByRole('button',{name:'Run segmentation'}).click();
  await page.getByRole('region',{name:'Depth vs allelic segments'}).getByRole('table',{name:'Breakpoint agreement'}).waitFor({timeout:900000});
  const runSeconds=Math.round((Date.now()-t0)/1000);
  const readTable=name=>page.getByRole('table',{name}).evaluate(t=>[...t.querySelectorAll('tbody tr')].map(r=>[...r.cells].map(c=>c.textContent.trim())));
  const chr8=await readTable('Breakpoint agreement');
  await page.getByRole('group',{name:'Comparison scope'}).getByRole('button',{name:'genome-wide'}).click();
  const genome=await readTable('Breakpoint agreement');
  const shapes=()=>page.evaluate(()=>{const s=document.querySelector('.allelic-view .js-plotly-plot').layout.shapes||[];const c={};for(const x of s)c[x.name||'other']=(c[x.name||'other']||0)+1;return c});
  const shapesAbs=await shapes();
  const [dl]=await Promise.all([page.waitForEvent('download'),page.getByRole('button',{name:'allelic_vs_depth.tsv'}).click()]);
  const tsv=fs.readFileSync(await dl.path(),'utf8');
  await page.getByLabel('Allelic signal',{exact:true}).selectOption('dbaf');
  await page.getByText(/Segments shown are from an earlier run/).waitFor();
  await page.waitForFunction(()=>{const e=document.querySelector('.allelic-view .js-plotly-plot');return !!(e&&e.layout&&e.layout.yaxis2&&e.layout.yaxis2.range&&e.layout.yaxis2.range[0]===-0.5)},null,{timeout:120000});
  const yRange=await page.evaluate(()=>document.querySelector('.allelic-view .js-plotly-plot').layout.yaxis2.range);
  await page.getByRole('group',{name:'Depth segments'}).getByRole('button',{name:'original'}).click();
  const out=path.join(__dirname,'..','validation','allelic');fs.mkdirSync(out,{recursive:true});
  await page.screenshot({path:path.join(out,'allelic.png'),fullPage:true});
  const result={panels,runSeconds,chr8,genome,shapesAbs,tsvHead:tsv.split('\n').slice(0,8),yRangeSigned:yRange,errors};
  fs.writeFileSync(path.join(out,'allelic_check_results.json'),JSON.stringify(result,null,1));
  assert.equal(panels,10,'5 clusters × (depth + allelic)');
  assert.ok(genome.length===5,'one agreement row per cluster');
  assert.ok(tsv.startsWith('# signal: abs_dbaf'),'export header');
  assert.deepEqual(yRange,[-0.5,0.5],'signed y range');
  assert.deepEqual(errors,[],'no page errors');
  console.log('PASS',JSON.stringify({runSeconds,genome}));
 } finally { await browser.close() }
})().catch(e=>{console.error(e);process.exit(1)});
