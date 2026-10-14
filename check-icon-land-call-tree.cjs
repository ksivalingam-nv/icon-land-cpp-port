// SPDX-License-Identifier: BSD-3-Clause
const fs = require('fs');
const path = require('path');
const {pathToFileURL} = require('url');
const { chromium } = require('playwright');
(async () => {
 console.log('Launching temporary browser');
 const launchOptions={headless:true};
 if(process.env.CALL_TREE_BROWSER_EXECUTABLE)launchOptions.executablePath=process.env.CALL_TREE_BROWSER_EXECUTABLE;
 const browser = await chromium.launch(launchOptions);
 const page = await browser.newPage({viewport:{width:1440,height:1100}});
 page.setDefaultTimeout(10000);
 const errors=[];page.on('pageerror',e=>{errors.push(e.message);console.log('PAGE ERROR',e.message)});
 console.log('Loading generated HTML');
 await page.goto(pathToFileURL(path.join(__dirname,'icon-land-call-tree.html')).href);
 await page.waitForSelector('.tree-row');
 await page.waitForFunction(()=>!!window.CALL_TREE_DEBUG);
 console.log('Page loaded');
 const diagnostics = await page.evaluate(() => {
  const {data}=window.CALL_TREE_DEBUG;
  const id=(name,mod)=>Object.values(data.nodes).find(n=>n.name.toLowerCase()===name.toLowerCase()&&n.module===mod&&n.kind!=='interface').id;
  const has=(from,to)=>data.edges[from]?.some(e=>e.to===to);
  const broken=[];for(const [from,edges] of Object.entries(data.edges))for(const e of edges)if(!data.nodes[e.to])broken.push(e.to);
  return {broken,stats:data.stats,taskEntries:data.tasks.length,distinctTasks:new Set(data.tasks.map(t=>t.name)).size,
   taskMappings:data.tasks.every(t=>data.nodes[t.integrate]&&data.nodes[t.aggregate]),
   entry:has(id('icon',''),id('jsbach_model','mo_jsbach_model')),
   timestep:has(id('jsbach_model','mo_jsbach_model'),id('run_one_timestep','mo_jsbach_model')),
   runTasks:has(id('run_one_timestep','mo_jsbach_model'),id('Run_tasks','mo_jsb_model_class')),
   doIt:has(id('Run_tasks','mo_jsb_model_class'),id('Do_it','mo_jsb_task_class')),
   dispatchTargets:(data.edges[Object.values(data.nodes).find(n=>n.name.toLowerCase()==='this%integrate'&&n.module==='mo_jsb_task_class').id]||[]).length};
 });
 if(diagnostics.broken.length||!diagnostics.taskMappings||!diagnostics.entry||!diagnostics.timestep||!diagnostics.runTasks||!diagnostics.doIt||diagnostics.dispatchTargets!==21)throw Error(JSON.stringify(diagnostics));
 await page.screenshot({path:path.join(__dirname,'icon-land-call-tree-desktop.png'),fullPage:true});
 console.log('Checking tasks');
 await page.getByRole('tab',{name:'Scientific tasks',exact:true}).click();
 const taskRows=await page.locator('#tree > ul > li').count();if(taskRows!==22)throw Error('Expected 22 task rows, got '+taskRows);
 await page.locator('#process-picker').selectOption('SSE_');
 if(await page.locator('#tree > ul > li').count()!==2)throw Error('SSE filter failed');
 await page.locator('#tree > ul > li').first().locator('summary').first().click();
 await page.waitForTimeout(100);
 await page.screenshot({path:path.join(__dirname,'icon-land-call-tree-tasks.png'),fullPage:true});
 await page.getByRole('tab',{name:'Kokkos bridge',exact:true}).click();
 await page.getByRole('tab',{name:'Routine explorer',exact:true}).click();
 console.log('Checking search');
 await page.getByRole('searchbox').fill('calc_vol_heat_capacity');
 await page.waitForTimeout(180);
 if(await page.locator('.result-button').count()<1)throw Error('Search failed');
 await page.locator('.result-button').first().click();
 if(!(await page.locator('#detail').innerText()).includes('calc_vol_heat_capacity'))throw Error('Selection failed');
 const outgoing=page.locator('#detail .edge-link').first();
 const target=await outgoing.getAttribute('href');
 await outgoing.click();
 if(await page.evaluate(()=>location.hash)!==target)throw Error('Routine link did not navigate');
 await page.reload();await page.waitForFunction(()=>!!window.CALL_TREE_DEBUG);
 if(await page.evaluate(()=>window.CALL_TREE_DEBUG.state.root)!==target.split('=')[1])throw Error('Deep link did not restore');
 const nativeContext=await browser.newContext({javaScriptEnabled:false,viewport:{width:1440,height:1100}});
 const native=await nativeContext.newPage();
 await native.goto(pathToFileURL(path.join(__dirname,'icon-land-call-tree.html')).href);
 if(await native.locator('#tree a.native-link').count()<15)throw Error('Native application tree is empty');
 const brokenNative=await native.evaluate(()=>[...document.querySelectorAll('a[href^="#"]')].filter(a=>!document.getElementById(a.hash.slice(1))).map(a=>a.hash));
 if(brokenNative.length)throw Error('Broken native targets: '+brokenNative.slice(0,8));
 await native.getByRole('tab',{name:'Scientific tasks',exact:true}).click();
 if(await native.evaluate(()=>location.hash)!=='#tasks')throw Error('Native tab link failed');
 const nativeTask=native.locator('#tasks li > details').first();await nativeTask.locator('summary').click();
 const nativeLink=nativeTask.locator('a').first();const nativeHref=await nativeLink.getAttribute('href');await nativeLink.click();
 if(await native.evaluate(()=>location.hash)!==nativeHref)throw Error('Native routine link failed');
 const nativeTarget=native.locator(nativeHref);await nativeTarget.locator('summary').click();
 if(!await nativeTarget.evaluate(d=>d.open))throw Error('Native source reference did not expand');
 await native.screenshot({path:path.join(__dirname,'icon-land-call-tree-no-javascript.png'),fullPage:false});
 await nativeContext.close();
 const blocked=await browser.newPage({viewport:{width:1440,height:1100}});
 await blocked.route('**/icon-land-call-tree-data.js',route=>route.abort());
 await blocked.goto(pathToFileURL(path.join(__dirname,'icon-land-call-tree.html')).href);
 if(await blocked.locator('#tree a.native-link').count()<15)throw Error('Blocked data file leaves an empty page');
 await blocked.getByRole('tab',{name:'Kokkos bridge',exact:true}).click();
 if(await blocked.evaluate(()=>location.hash)!=='#bridge')throw Error('Blocked data fallback navigation failed');
 await blocked.close();
 await page.locator('#utilities').check();
 await page.locator('#overview').click();
 await page.getByRole('tab',{name:'Application',exact:true}).click();
 await page.setViewportSize({width:390,height:844});
 await page.screenshot({path:path.join(__dirname,'icon-land-call-tree-mobile.png'),fullPage:true});
 const overflow=await page.evaluate(()=>document.documentElement.scrollWidth>window.innerWidth);
 if(overflow)throw Error('Mobile horizontal overflow');
 await page.emulateMedia({colorScheme:'dark'});
 await page.setViewportSize({width:1440,height:1100});
 await page.screenshot({path:path.join(__dirname,'icon-land-call-tree-dark.png'),fullPage:true});
 if(errors.length)throw Error('Browser errors: '+errors.join('; '));
 console.log(JSON.stringify({...diagnostics,browserErrors:errors,taskRows,mobileOverflow:overflow,nativeLinks:'passed with JavaScript disabled',blockedData:'usable fallback',deepLinks:'passed'},null,2));
 await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
