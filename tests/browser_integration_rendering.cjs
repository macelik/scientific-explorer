const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:true,executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 await page.addInitScript(()=>{window.plotRenderCount=0;window.longTasks=[];new PerformanceObserver(l=>window.longTasks.push(...l.getEntries().map(e=>e.duration))).observe({entryTypes:['longtask']});});
 await page.goto('http://127.0.0.1:8766');await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});
 const offscreen=await page.locator('#cluster-tracks .js-plotly-plot').count();
 assert.equal(offscreen,0,'Do not initialize offscreen chromosome Plotly charts during CN-panel load');
 const tracks=page.locator('#cluster-tracks');await tracks.scrollIntoViewIfNeeded();
 const plot=tracks.locator('.js-plotly-plot');await plot.waitFor();
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?._fullLayout);
 const axes=await plot.evaluate(el=>({traces:el.data.length,shapes:el.layout.shapes.length,range:el._fullLayout.xaxis.range}));
 assert(axes.traces>=5);assert(axes.shapes>0);
 await page.getByRole('checkbox',{name:'allelic |dBAF| pseudobulk (right axis)',exact:true}).check();
 await page.getByLabel(/^track level/).selectOption('raw');
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?.data.filter(t=>t.name?.includes('allelic')).length===5);
 const diagnostics=page.locator('.panel').filter({has:page.getByRole('heading',{name:/Small-segment diagnostics across/})});
 const slots=diagnostics.locator('.deferred-plot');assert.equal(await slots.count(),7);
 for(let i=0;i<7;i++){await slots.nth(i).scrollIntoViewIfNeeded();await slots.nth(i).locator('.js-plotly-plot').waitFor();}
 assert.equal(await diagnostics.locator('.js-plotly-plot').count(),7);
 await page.getByRole('button',{name:'Cohort karyogram',exact:true}).click();
 await page.getByText('300 matching cells',{exact:true}).waitFor();
 const out=path.resolve(__dirname,'../validation/integration-freeze');fs.mkdirSync(out,{recursive:true});
 fs.writeFileSync(path.join(out,'rendering-results.json'),JSON.stringify({offscreen_charts_on_load:offscreen,axes,diagnostics_loaded:7,controls:true,navigation:true,errors},null,2));
 assert.deepEqual(errors,[]);console.log('PASS deferred segmentation plots, tracks, allelic controls, all diagnostics and navigation');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
