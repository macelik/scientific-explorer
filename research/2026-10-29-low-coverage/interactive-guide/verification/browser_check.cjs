const { chromium } = require(process.env.PLAYWRIGHT_MODULE || '/home/mami/.npm/_npx/9833c18b2d85bc59/node_modules/playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const { pathToFileURL } = require('node:url');
const root = path.resolve(__dirname, '..');
const checks = [], errors = [], network = [];
const check = (name, passed) => { assert.ok(passed, name); checks.push(name); };
async function capture(page, name) {
  await page.evaluate(() => window.scrollTo(0, 0));
  await page.screenshot({ path: path.join(__dirname, name), fullPage: true, animations: 'disabled' });
}

(async () => {
  const browser = await chromium.launch({
    executablePath: process.env.CHROMIUM_PATH || '/home/mami/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome',
    headless: true,
  });
  const context = await browser.newContext({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
  const page = await context.newPage();
  page.on('pageerror', e => errors.push(String(e)));
  page.on('request', r => { if (/^https?:/.test(r.url())) network.push(r.url()); });
  const url = pathToFileURL(path.join(root, 'index.html')).href;
  await page.goto(url);
  await page.waitForSelector('h1');
  const data = await page.locator('#guide-data').textContent().then(JSON.parse);
  check('36 methods loaded offline', data.methods.length === 36 && await page.locator('.method-button').count() === 36);
  await capture(page, 'desktop-method.png');

  for (const m of data.methods) {
    console.log(`Checking ${m.id}: ${m.name}`);
    await page.locator(`.method-button[data-method="${m.id}"]`).click();
    check(`${m.id}: heading and steps`, (await page.locator('h1').textContent()) === m.name && await page.locator('.stage-node').count() === m.steps.length);
    for (let i = 0; i < m.steps.length; i++) {
      await page.locator(`.stage-node[data-step="${i}"]`).click();
      assert.equal(await page.locator('.step-focus p').textContent(), m.steps[i].text);
    }
    check(`${m.id}: step retains keyboard focus`, await page.locator('.stage-node.active').evaluate(e => e === document.activeElement));
    await page.locator('[data-panel="example"]').click();
    check(`${m.id}: worked example`, (await page.locator('#main').textContent()).includes(m.example.interpretation));
    await page.locator('[data-panel="adapt"]').click();
    check(`${m.id}: workflow adaptation`, (await page.locator('#main').textContent()).includes(m.adaptation.first_experiment));
  }

  await page.locator('#search').fill('Numbat');
  check('Search filters methods', (await page.locator('.method-button').count()) > 0 && (await page.locator('.method-button').count()) < 36);
  await page.locator('#search').fill('no-such-method-xyz');
  check('Empty search is recoverable', await page.locator('[data-action="clear-search"]').isVisible());
  await page.locator('[data-action="clear-search"]').click();
  for (const [cat, n] of Object.entries({ A: 11, B: 5, C: 6, D: 6, E: 3, F: 5 })) {
    await page.locator('#category').selectOption(cat);
    check(`Category ${cat}`, await page.locator('.method-button').count() === n);
  }
  await page.locator('#category').selectOption('all');
  await page.locator('[data-view="compare"]').click();
  await page.locator('[data-action="suggest-compare"]').click();
  check('Three-method comparison', await page.locator('.compare-column').count() === 3);
  await capture(page, 'desktop-compare.png');
  await page.locator('.method-button[data-method="A2"]').click();
  await page.locator('[data-pin="A2"]').click();
  check('Comparison limit enforced', await page.locator('#compare-count').textContent() === '3' && (await page.locator('#toast').textContent()).includes('up to three'));
  await page.locator('[data-view="compare"]').click();
  await page.locator('[data-pin="A1"]').click();
  check('Comparison removal', await page.locator('.compare-column').count() === 2);
  await page.locator('[data-action="clear-compare"]').click();
  check('Comparison clear', await page.locator('#compare-count').textContent() === '0');

  await page.locator('[data-view="examples"]').click();
  for (const lab of ['depth', 'reference', 'alleles', 'zeros']) {
    await page.locator(`.lab-tabs [data-lab="${lab}"]`).click();
    check(`Lab ${lab} chart`, await page.locator('#lab-chart svg').count() === 1);
    if (lab === 'depth') {
      await page.locator('[data-preset="copy"]').click();
      const a = await page.locator('#lab-chart path').getAttribute('d');
      await page.locator('[data-preset="access"]').click();
      check('Copy/accessibility confounding gives identical curve', a === await page.locator('#lab-chart path').getAttribute('d'));
      await page.locator('#access').focus();
      await page.keyboard.press('ArrowRight');
      check('Keyboard slider changes curve and retains focus', a !== await page.locator('#lab-chart path').getAttribute('d') && await page.locator('#access').evaluate(e => e === document.activeElement));
    }
    if (lab === 'reference') {
      await page.locator('[data-preset="copy"]').click();
      const a = await page.locator('#lab-readouts').textContent();
      await page.locator('[data-preset="reference"]').click();
      check('Reference choice changes ratio', a !== await page.locator('#lab-readouts').textContent());
    }
    if (lab === 'alleles') {
      await page.locator('[data-preset="no-alleles"]').click();
      check('No alleles remains undefined and prior-only', /Undefined.*Prior only/s.test(await page.locator('#lab-readouts').textContent()));
      await page.locator('[data-preset="many-alleles"]').click();
      await capture(page, 'desktop-alleles.png');
    }
    if (lab === 'zeros') {
      for (const model of ['poisson', 'nb', 'zinb']) {
        await page.locator('#count-model').selectOption(model);
        check(`Count model ${model}: controls`, await page.locator('#r').isDisabled() === (model === 'poisson') && await page.locator('#pi').isDisabled() === (model !== 'zinb'));
      }
    }
  }

  const fixtures = JSON.parse(fs.readFileSync(path.join(__dirname, 'math_reference.json'), 'utf8'));
  const maths = await page.evaluate(fixtures => ({
    counts: fixtures.counts.map(x => GuideMath.countDistribution(x.model, x.mu, x.r, x.pi)),
    alleles: fixtures.alleles.map(x => GuideMath.allelePosterior(x.n, x.hb)),
  }), fixtures);
  fixtures.counts.forEach((x, i) => check(`Distribution vs SciPy: ${x.model}/${x.mu}/${x.r}/${x.pi}`, x.probabilities.every((p, k) => Math.abs(p - maths.counts[i][k]) < 1e-11)));
  fixtures.alleles.forEach((x, i) => check(`Beta interval vs SciPy: ${x.n}/${x.hb}`, Math.abs(x.mean - maths.alleles[i].mean) < 1e-12 && Math.abs(x.low - maths.alleles[i].low) < 0.0004 && Math.abs(x.high - maths.alleles[i].high) < 0.0004));

  await page.locator('[data-view="pilot"]').click();
  for (const row of data.pilot) {
    await page.locator('#pilot-region').selectOption(row.region);
    await page.locator('#pilot-stratum').selectOption(`${row.stratum}:${row.level}`);
    const summary = await page.locator('.pilot-summary strong').allTextContents();
    const bars = await page.locator('.bar-row > span:last-child').allTextContents();
    check(`Pilot ${row.region}/${row.stratum}/${row.level}`, Number(summary[0]) === row.n_cells && Number(summary[1]) === row.n_bins && bars[0] === `${row.raw_whole_zero_cells} / ${row.n_cells}` && bars[1] === `${row.at_whole_zero_cells} / ${row.n_cells}`);
  }
  await page.locator('#pilot-stratum').selectOption('all:all');
  await capture(page, 'desktop-pilot.png');
  await page.locator('[data-view="workflow"]').click();
  for (const stage of ['inputs', 'reference', 'segment', 'alleles', 'interpret', 'validate']) {
    await page.locator(`.workflow-step[data-workflow="${stage}"]`).click();
    check(`Workflow ${stage}`, await page.locator(`.workflow-step[data-workflow="${stage}"]`).getAttribute('aria-pressed') === 'true');
  }
  await page.locator('.glossary-open').click();
  check('Glossary modal opens', await page.locator('#glossary').isVisible());
  await page.keyboard.press('Escape');
  check('Glossary closes with Escape', !await page.locator('#glossary').isVisible());

  await page.goto(`${url}#view=methods&id=A4&panel=adapt`);
  check('Deep link restores method and panel', (await page.locator('h1').textContent()).includes('Numbat') && await page.locator('[data-panel="adapt"]').getAttribute('aria-pressed') === 'true');
  const [download] = await Promise.all([page.waitForEvent('download'), page.locator('[data-action="export"]').click()]);
  await download.saveAs(path.join(__dirname, 'example-export.md'));
  const exported = fs.readFileSync(path.join(__dirname, 'example-export.md'), 'utf8');
  check('Notes download contains mechanics, proposal and citations', ['## How it works', '## Proposed adaptation', 'https://'].every(x => exported.includes(x)));
  await page.emulateMedia({ media: 'print' });
  check('Print hides navigation', !await page.locator('.masthead').isVisible() && !await page.locator('.sidebar').isVisible());
  await page.emulateMedia({ media: 'screen' });

  const widths = [1440, 1024, 768, 390, 320];
  for (const width of widths) {
    await page.setViewportSize({ width, height: 900 });
    for (const view of ['methods', 'examples', 'pilot', 'workflow', 'compare']) {
      await page.locator(`[data-view="${view}"]`).click();
      const size = await page.evaluate(() => ({ scroll: document.documentElement.scrollWidth, client: document.documentElement.clientWidth }));
      check(`No page overflow ${width}px/${view}`, size.scroll <= size.client + 1);
    }
  }
  await page.setViewportSize({ width: 390, height: 844 });
  await page.locator('.method-button[data-method="A4"]').click();
  await capture(page, 'mobile-method.png');
  await page.locator('[data-view="examples"]').click();
  await page.locator('.lab-tabs [data-lab="zeros"]').click();
  await capture(page, 'mobile-example.png');
  check('No runtime errors', errors.length === 0);
  check('No external network requests', network.length === 0);
  fs.writeFileSync(path.join(__dirname, 'browser_report.json'), JSON.stringify({ status: 'passed', tested_html_sha256: require('node:crypto').createHash('sha256').update(fs.readFileSync(path.join(root, 'index.html'))).digest('hex'), browser: await browser.version(), check_count: checks.length, checks, errors, external_requests: network, viewports: widths, scope: 'Local file, real Chromium interactions, scientific demo math compared with SciPy; no CN caller was executed.' }, null, 2) + '\n');
  console.log(JSON.stringify({ status: 'passed', checks: checks.length, errors, external_requests: network }, null, 2));
  await browser.close();
})().catch(error => {
  fs.writeFileSync(path.join(__dirname, 'browser_report.json'), JSON.stringify({ status: 'failed', passed_checks: checks, errors, failure: String(error) }, null, 2) + '\n');
  console.error(error);
  process.exit(1);
});
