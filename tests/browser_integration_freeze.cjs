const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 const watchdog=setTimeout(()=>{console.error('FAIL browser responsiveness deadline');browser.close().finally(()=>process.exit(1))},90000);
 const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[],samples=[];
 page.on('pageerror',e=>errors.push(String(e)));
 await page.addInitScript(()=>{window.probe={ticks:0,long:[]};setInterval(()=>window.probe.ticks++,100);new PerformanceObserver(l=>window.probe.long.push(...l.getEntries().map(x=>({start:x.startTime,ms:x.duration})))).observe({entryTypes:['longtask']});});
 const cdp=await page.context().newCDPSession(page);await cdp.send('Profiler.enable');await cdp.send('Profiler.start');
 await page.goto('http://127.0.0.1:8766');await page.getByText('300 matching cells',{exact:true}).waitFor();
 const start=Date.now();await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});
 console.log('CN panels ready in',Date.now()-start,'ms');
 await page.locator('#cluster-tracks').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?._fullLayout);
 console.log('Tracks ready in',Date.now()-start,'ms');
 for(const width of [1280,1000,800,1920,1280]){
  await page.setViewportSize({width,height:900});
  await new Promise(r=>setTimeout(r,1500));
  const s=await page.evaluate(()=>({width:innerWidth,documentWidth:document.documentElement.scrollWidth,ticks:window.probe.ticks,long:window.probe.long.slice(-5),plotWidths:[...document.querySelectorAll('.js-plotly-plot')].map(e=>({width:e.clientWidth,svg:e._fullLayout?.width,height:e.clientHeight}))}));samples.push(s);console.log(JSON.stringify(s));
 }
 await page.getByRole('checkbox',{name:'allelic |dBAF| pseudobulk (right axis)',exact:true}).check();
 await page.getByLabel(/^track level/).selectOption('raw');
 await page.getByRole('button',{name:'Cohort karyogram',exact:true}).click();
 await page.getByText('300 matching cells',{exact:true}).waitFor();
 const out=path.resolve(__dirname,'../validation/integration-freeze');fs.mkdirSync(out,{recursive:true});
 fs.writeFileSync(path.join(out,'profile.json'),JSON.stringify((await cdp.send('Profiler.stop')).profile));
 fs.writeFileSync(path.join(out,'browser.json'),JSON.stringify({samples,errors},null,2));
 if(errors.length)throw Error(errors.join('\n'));
 clearTimeout(watchdog);await browser.close();console.log('PASS integration load, resize, controls and navigation');
})().catch(e=>{console.error(e);process.exit(1)});
