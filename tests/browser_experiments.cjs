const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs');
const path = require('path');
const out = path.resolve(__dirname, '../validation');

(async () => {
  const browser = await chromium.launch({ headless: true, executablePath: process.env.CHROMIUM_PATH });
  const page = await browser.newPage({ viewport: { width: 1500, height: 1080 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  const timings = {};
  const start = Date.now();
  for (let attempt=0; attempt<30; attempt++) {
    try { const response=await page.request.get('http://127.0.0.1:8766/api/health'); if(response.ok()) break; } catch {}
    await new Promise(r=>setTimeout(r,1000));
  }
  page.on('console', msg => { if (msg.type()==='error') console.error('BROWSER:', msg.text()); });
  const shot = async name => { await page.evaluate(()=>window.scrollTo(0,0)); await page.screenshot({path:path.join(out,name),fullPage:true}); };
  await page.goto('http://127.0.0.1:8766');
  await page.getByText('300 matching cells', { exact: true }).waitFor({ timeout: 30000 });
  timings.boot_ms = Date.now() - start;
  await shot('cohort.png');
  await page.getByRole('button', { name: 'TTAGGCTAGGCCGGAA-1 / chr9', exact: true }).click();
  await page.getByText('Signal → candidates → recursion', { exact: true }).waitFor();
  await shot('explorer.png');
  await page.getByRole('button', { name: 'X experiments', exact: true }).click();
  await page.getByRole('button', { name: 'Run experiment', exact: true }).waitFor();
  const run = async (name) => {
    const begin = Date.now();
    await page.getByRole('button', { name: 'Run experiment', exact: true }).click();
    await page.getByRole('button', { name: 'Run experiment', exact: true }).waitFor({ timeout: 90000 });
    await page.getByRole('heading', { name: 'Original → experiment', exact: true }).waitFor();
    timings[name] = Date.now() - begin;
  };
  await run('multiply_ms');
  await page.getByRole('button', { name: 'Export experiment JSON', exact: true }).click();
  await shot('experiment-multiply.png');
  await page.getByLabel('X transformation', { exact: true }).selectOption('simulate');
  await page.getByLabel('Simulation class', { exact: true }).selectOption('loss');
  await run('simulate_ms');
  await page.getByText('Both children and exact split tests', { exact: true }).click();
  await shot('experiment-simulate.png');
  await page.getByLabel('X transformation', { exact: true }).selectOption('duplicate');
  await page.getByLabel('Extension direction', { exact: true }).selectOption('left');
  await page.getByLabel('Extension bins', { exact: true }).fill('50');
  await run('duplicate_ms');
  await page.getByLabel('X transformation', { exact: true }).selectOption('extend');
  await page.getByLabel('Extension direction', { exact: true }).selectOption('right');
  await page.getByLabel('Extension bins', { exact: true }).fill('50');
  await run('extend_ms');
  await page.getByRole('button', { name: 'Run width series', exact: true }).click();
  await page.getByRole('button', { name: 'Stop after current run', exact: true }).waitFor();
  await page.getByRole('button', { name: 'Stop after current run', exact: true }).waitFor({ state: 'hidden', timeout: 120000 });
  await page.getByRole('button', { name: 'Inspect run 9', exact: true }).waitFor();
  await shot('width-series.png');
  // Editing marks frozen results stale; it does not pretend old AD is recomputed.
  await page.getByLabel('Extension bins', { exact: true }).fill('51');
  await page.getByText(/Settings changed\. These are the previous run/).waitFor();
  await page.getByRole('button', { name: 'Explore output', exact: true }).click();
  await page.getByRole('button', { name: 'X experiments', exact: true }).click();
  if (await page.getByLabel('Extension bins', { exact: true }).inputValue() !== '51') throw new Error('Draft lost after navigation');
  // Persist the experiment history through the shared session serializer.
  const downloadPromise = page.waitForEvent('download');
  await page.getByRole('button', { name: 'download JSON', exact: true }).click();
  const download = await downloadPromise;
  const sessionPath = path.join(out, 'browser-session.json');
  await download.saveAs(sessionPath);
  const session = JSON.parse(fs.readFileSync(sessionPath, 'utf8'));
  const exp = session.state.experiments['TTAGGCTAGGCCGGAA-1|chr9'];
  if (exp.history.length !== 9 || exp.draft.length !== 51) throw new Error('Experiment history missing from session');
  await page.reload();
  try { await page.getByText('300 matching cells', { exact: true }).waitFor(); }
  catch (e) {
    await shot('reload-failure.png');
    fs.writeFileSync(path.join(out,'reload-failure.txt'), `${page.url()}\n${await page.locator('body').innerText()}\n${errors.join('\n')}`);
    throw e;
  }
  await page.locator('input[type=file]').setInputFiles(sessionPath);
  await page.getByRole('heading', { name: 'X-value experiments', exact: true }).waitFor();
  await page.getByRole('button', { name: 'Inspect run 9', exact: true }).waitFor();
  if (await page.getByLabel('Extension bins', { exact: true }).inputValue() !== '51') throw new Error('Draft not restored from session');
  await page.setViewportSize({ width: 1050, height: 900 });
  await page.waitForFunction(() => document.documentElement.scrollWidth <= window.innerWidth + 2);
  await page.waitForFunction(() => [...document.querySelectorAll('.experiment-preview > .js-plotly-plot, .experiment-results > .js-plotly-plot')].every(el => el.getBoundingClientRect().height >= el._fullLayout.height));
  await shot('experiment-narrow.png');
  if (errors.length) throw new Error(errors.join('\n'));
  fs.writeFileSync(path.join(out, 'browser-results.json'), JSON.stringify({ timings, errors, history_runs: 9, session_restored: true }, null, 2));
  console.log(JSON.stringify({ timings, errors, history_runs: 9, session_restored: true }));
  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
