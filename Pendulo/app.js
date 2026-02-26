// Péndulo simple - Simulador (canvas + RK4 + arrastre + plots + CSV + T medido/ciclos + tabla)

const simCanvas = document.getElementById("simCanvas");
const ctx = simCanvas.getContext("2d");

const thetaPlot = document.getElementById("thetaPlot");
const tctx = thetaPlot.getContext("2d");
const omegaPlot = document.getElementById("omegaPlot");
const octx = omegaPlot.getContext("2d");

// UI
const Ls = document.getElementById("L");
const ms = document.getElementById("m");
const th0s = document.getElementById("th0");
const gs = document.getElementById("g");
const betas = document.getElementById("beta");
const nonlinearToggle = document.getElementById("nonlinearToggle");
const dampingToggle = document.getElementById("dampingToggle");

const LVal = document.getElementById("LVal");
const mVal = document.getElementById("mVal");
const th0Val = document.getElementById("th0Val");
const gVal = document.getElementById("gVal");
const betaVal = document.getElementById("betaVal");

const Ttheory = document.getElementById("Ttheory");
const ftheory = document.getElementById("ftheory");
const thetaNow = document.getElementById("thetaNow");
const omegaNow = document.getElementById("omegaNow");
const Enow = document.getElementById("Enow");

const TmeasuredEl = document.getElementById("Tmeasured");
const cyclesCountEl = document.getElementById("cyclesCount");
const TavgEl = document.getElementById("Tavg");

const dataBody = document.getElementById("dataBody");

const btnStart = document.getElementById("btnStart");
const btnReset = document.getElementById("btnReset");
const btnResetMeasure = document.getElementById("btnResetMeasure");
const btnCSV = document.getElementById("btnCSV");

// Sim state
let running = false;
let dragging = false;

let t = 0;
let theta = deg2rad(+th0s.value); // rad
let omega = 0;                    // rad/s

// data for plots + CSV
const maxSeconds = 20;
const dt = 1 / 120; // 120 Hz
let series = [];    // {t, theta, omega, alpha}

// --- Medición experimental de período ---
let lastCrossTime = null;     // instante del último cruce válido
let periodMeasured = null;    // último T medido
let periods = [];             // lista de periodos recientes (máx 10)
let cyclesCount = 0;          // número de ciclos contados

function resetMeasurement() {
  lastCrossTime = null;
  periodMeasured = null;
  periods = [];
  cyclesCount = 0;
}

// Resize canvases to match CSS size
function fitCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = Math.max(2, Math.floor(rect.width * dpr));
  canvas.height = Math.max(2, Math.floor(rect.height * dpr));
  const c = canvas.getContext("2d");
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { w: rect.width, h: rect.height };
}

function syncAllCanvasSizes() {
  fitCanvas(simCanvas);
  fitCanvas(thetaPlot);
  fitCanvas(omegaPlot);
}

window.addEventListener("resize", () => {
  syncAllCanvasSizes();
  draw();
  drawPlots();
  renderTable();
});

// Physics helpers
function deg2rad(d) { return (d * Math.PI) / 180; }
function rad2deg(r) { return (r * 180) / Math.PI; }
function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }

function params() {
  const L = +Ls.value;
  const m = +ms.value;
  const g = +gs.value;
  const beta = dampingToggle.checked ? +betas.value : 0;
  const nonlinear = nonlinearToggle.checked;
  return { L, m, g, beta, nonlinear };
}

// ODE: theta' = omega
//      omega' = -2 beta omega - (g/L) * sin(theta)   (nonlinear)
//      omega' = -2 beta omega - (g/L) * theta        (small-angle)
function alpha(th, om) {
  const { L, g, beta, nonlinear } = params();
  const restoring = nonlinear ? Math.sin(th) : th;
  return -2 * beta * om - (g / L) * restoring;
}

// RK4 step
function stepRK4() {
  const prevTheta = theta;
  const prevOmega = omega;

  const th = theta, om = omega;

  const k1_th = om;
  const k1_om = alpha(th, om);

  const k2_th = om + 0.5 * dt * k1_om;
  const k2_om = alpha(th + 0.5 * dt * k1_th, om + 0.5 * dt * k1_om);

  const k3_th = om + 0.5 * dt * k2_om;
  const k3_om = alpha(th + 0.5 * dt * k2_th, om + 0.5 * dt * k2_om);

  const k4_th = om + dt * k3_om;
  const k4_om = alpha(th + dt * k3_th, om + dt * k3_om);

  theta = th + (dt / 6) * (k1_th + 2 * k2_th + 2 * k3_th + k4_th);
  omega = om + (dt / 6) * (k1_om + 2 * k2_om + 2 * k3_om + k4_om);

  t += dt;

  const a = alpha(theta, omega);
  series.push({ t, theta, omega, alpha: a });

  // keep only last maxSeconds
  const tMin = t - maxSeconds;
  while (series.length && series[0].t < tMin) series.shift();

  // --- Detectar cruce por equilibrio: de theta<0 a theta>=0 con omega>0 (1 vez por período) ---
  if (prevTheta < 0 && theta >= 0 && omega > 0) {
    if (lastCrossTime !== null) {
      const T = t - lastCrossTime;
      // descarta valores raros por “rebotes” numéricos
      if (T > 0.1 && T < 20) {
        periodMeasured = T;
        periods.push(T);
        if (periods.length > 10) periods.shift();
        cyclesCount += 1;
      }
    }
    lastCrossTime = t;
  }
}

// UI updates
function updateLabels() {
  const { L, m, g, beta, nonlinear } = params();
  LVal.textContent = L.toFixed(2);
  mVal.textContent = m.toFixed(2);
  th0Val.textContent = (+th0s.value).toFixed(0);
  gVal.textContent = g.toFixed(2);
  betaVal.textContent = beta.toFixed(2);

  const T = 2 * Math.PI * Math.sqrt(L / g);
  const f = 1 / T;
  Ttheory.textContent = `${T.toFixed(3)} s`;
  ftheory.textContent = `${f.toFixed(3)} Hz`;

  thetaNow.textContent = `${rad2deg(theta).toFixed(2)} °`;
  omegaNow.textContent = `${omega.toFixed(3)} rad/s`;

  // Energy
  const U = nonlinear ? m * g * L * (1 - Math.cos(theta)) : 0.5 * m * g * L * theta * theta;
  const K = 0.5 * m * (L * omega) * (L * omega);
  const E = K + U;
  Enow.textContent = `${E.toFixed(4)} J`;

  if (periodMeasured === null) TmeasuredEl.textContent = "—";
  else TmeasuredEl.textContent = `${periodMeasured.toFixed(3)} s`;

  cyclesCountEl.textContent = `${cyclesCount}`;

  if (!periods.length) {
    TavgEl.textContent = "—";
  } else {
    const avg = periods.reduce((s, x) => s + x, 0) / periods.length;
    TavgEl.textContent = `${avg.toFixed(3)} s`;
  }
}

function resetState(useSliderAngle = true) {
  t = 0;
  theta = useSliderAngle ? deg2rad(+th0s.value) : theta;
  omega = 0;

  resetMeasurement();

  series = [{ t: 0, theta, omega, alpha: alpha(theta, omega) }];
  draw();
  drawPlots();
  renderTable();
  updateLabels();
}

// Drawing pendulum
function draw() {
  const rect = simCanvas.getBoundingClientRect();
  const w = rect.width, h = rect.height;

  ctx.clearRect(0, 0, w, h);

  // coordinate system
  const pivot = { x: w * 0.5, y: 60 };
  const { L } = params();

  // scale L meters to pixels
  const pxPerMeter = Math.min(220, (h - 140) / 2.2);
  const rod = L * pxPerMeter;

  // bob position
  const x = pivot.x + rod * Math.sin(theta);
  const y = pivot.y + rod * Math.cos(theta);

  // background grid
  ctx.save();
  ctx.globalAlpha = 0.18;
  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 1;
  const step = 40;
  for (let gx = 0; gx <= w; gx += step) { ctx.beginPath(); ctx.moveTo(gx, 0); ctx.lineTo(gx, h); ctx.stroke(); }
  for (let gy = 0; gy <= h; gy += step) { ctx.beginPath(); ctx.moveTo(0, gy); ctx.lineTo(w, gy); ctx.stroke(); }
  ctx.restore();

  // ceiling bar
  ctx.save();
  ctx.strokeStyle = "rgba(255,255,255,.35)";
  ctx.lineWidth = 6;
  ctx.beginPath();
  ctx.moveTo(pivot.x - 90, pivot.y - 22);
  ctx.lineTo(pivot.x + 90, pivot.y - 22);
  ctx.stroke();
  ctx.restore();

  // rod
  ctx.save();
  ctx.strokeStyle = "rgba(255,255,255,.75)";
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.moveTo(pivot.x, pivot.y);
  ctx.lineTo(x, y);
  ctx.stroke();
  ctx.restore();

  // pivot point
  ctx.save();
  ctx.fillStyle = "#22c55e";
  ctx.beginPath();
  ctx.arc(pivot.x, pivot.y, 6, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();

  // bob
  const { m } = params();
  const r = clamp(12 + 10 * Math.sqrt(m), 12, 26);

  ctx.save();
  ctx.fillStyle = "rgba(34,197,94,.25)";
  ctx.strokeStyle = "rgba(34,197,94,.85)";
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, Math.PI * 2);
  ctx.fill();
  ctx.stroke();
  ctx.restore();

  // angle indicator + HUD
  ctx.save();
  ctx.strokeStyle = "rgba(255,255,255,.35)";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(pivot.x, pivot.y);
  ctx.lineTo(pivot.x, pivot.y + rod);
  ctx.stroke();

  ctx.fillStyle = "rgba(255,255,255,.7)";
  ctx.font = "12px system-ui, sans-serif";
  ctx.fillText(`θ = ${rad2deg(theta).toFixed(1)}°`, 12, 22);
  ctx.fillText(`t = ${t.toFixed(2)} s`, 12, 40);
  ctx.restore();

  // store for dragging
  lastGeom = { pivot, x, y, r };
}

let lastGeom = null;

// Plots (simple canvas plotting)
function drawPlots() {
  drawTimeSeries(thetaPlot, tctx, series.map(p => ({ x: p.t, y: rad2deg(p.theta) })), "θ (°)");
  drawTimeSeries(omegaPlot, octx, series.map(p => ({ x: p.t, y: p.omega })), "ω (rad/s)");
}

function drawTimeSeries(canvas, c, pts, yLabel) {
  const rect = canvas.getBoundingClientRect();
  const w = rect.width, h = rect.height;

  c.clearRect(0, 0, w, h);

  // frame
  c.strokeStyle = "rgba(255,255,255,.12)";
  c.lineWidth = 1;
  c.strokeRect(10, 10, w - 20, h - 20);

  if (pts.length < 2) return;

  const xmin = pts[0].x;
  const xmax = pts[pts.length - 1].x;

  let ymin = Infinity, ymax = -Infinity;
  for (const p of pts) { ymin = Math.min(ymin, p.y); ymax = Math.max(ymax, p.y); }
  if (Math.abs(ymax - ymin) < 1e-9) { ymax += 1; ymin -= 1; }

  const pad = 18;
  const X = x => 10 + pad + ((x - xmin) / (xmax - xmin)) * (w - 20 - 2 * pad);
  const Y = y => 10 + pad + (1 - (y - ymin) / (ymax - ymin)) * (h - 20 - 2 * pad);

  // grid
  c.save();
  c.globalAlpha = 0.35;
  c.strokeStyle = "rgba(255,255,255,.10)";
  for (let i = 1; i <= 4; i++) {
    const yy = 10 + pad + (i / 5) * (h - 20 - 2 * pad);
    c.beginPath(); c.moveTo(10 + pad, yy); c.lineTo(w - 10 - pad, yy); c.stroke();
  }
  c.restore();

  // plot
  c.save();
  c.strokeStyle = "rgba(34,197,94,.85)";
  c.lineWidth = 2;
  c.beginPath();
  c.moveTo(X(pts[0].x), Y(pts[0].y));
  for (let i = 1; i < pts.length; i++) c.lineTo(X(pts[i].x), Y(pts[i].y));
  c.stroke();
  c.restore();

  // labels
  c.fillStyle = "rgba(255,255,255,.70)";
  c.font = "12px system-ui, sans-serif";
  c.fillText(yLabel, 14, 20);
  c.fillText(`${xmin.toFixed(1)}s`, 14, h - 12);
  c.fillText(`${xmax.toFixed(1)}s`, w - 62, h - 12);
}

// Tabla (últimas 10 muestras)
function renderTable() {
  const last = series.slice(-10);
  let html = "";
  for (const p of last) {
    html += `<tr>
      <td>${p.t.toFixed(2)}</td>
      <td>${rad2deg(p.theta).toFixed(2)}</td>
      <td>${p.omega.toFixed(3)}</td>
      <td>${periodMeasured ? periodMeasured.toFixed(3) : "—"}</td>
      <td>${cyclesCount}</td>
    </tr>`;
  }
  dataBody.innerHTML = html;
}

// Interaction: drag bob to set theta0
function pointerPos(evt) {
  const rect = simCanvas.getBoundingClientRect();
  const x = (evt.touches ? evt.touches[0].clientX : evt.clientX) - rect.left;
  const y = (evt.touches ? evt.touches[0].clientY : evt.clientY) - rect.top;
  return { x, y };
}

function hitBob(p) {
  if (!lastGeom) return false;
  const dx = p.x - lastGeom.x;
  const dy = p.y - lastGeom.y;
  return Math.hypot(dx, dy) <= lastGeom.r + 8;
}

function setThetaFromPointer(p) {
  if (!lastGeom) return;
  const { pivot } = lastGeom;
  const dx = p.x - pivot.x;
  const dy = p.y - pivot.y;
  theta = Math.atan2(dx, dy);
  theta = clamp(theta, deg2rad(-85), deg2rad(85));
  omega = 0;

  // actualiza slider θ0 con magnitud del ángulo
  const d = Math.abs(rad2deg(theta));
  th0s.value = clamp(Math.round(d), +th0s.min, +th0s.max);
}

simCanvas.addEventListener("mousedown", (e) => {
  const p = pointerPos(e);
  if (hitBob(p)) {
    dragging = true;
    running = false;
    btnStart.textContent = "Iniciar";
    setThetaFromPointer(p);
    resetState(false);
  }
});

simCanvas.addEventListener("mousemove", (e) => {
  if (!dragging) return;
  setThetaFromPointer(pointerPos(e));
  draw();
  updateLabels();
});

window.addEventListener("mouseup", () => { dragging = false; });

simCanvas.addEventListener("touchstart", (e) => {
  const p = pointerPos(e);
  if (hitBob(p)) {
    dragging = true;
    running = false;
    btnStart.textContent = "Iniciar";
    setThetaFromPointer(p);
    resetState(false);
  }
}, { passive: true });

simCanvas.addEventListener("touchmove", (e) => {
  if (!dragging) return;
  setThetaFromPointer(pointerPos(e));
  draw();
  updateLabels();
}, { passive: true });

window.addEventListener("touchend", () => { dragging = false; });

// Buttons
btnStart.addEventListener("click", () => {
  running = !running;
  btnStart.textContent = running ? "Pausar" : "Iniciar";
});

btnReset.addEventListener("click", () => {
  running = false;
  btnStart.textContent = "Iniciar";
  resetState(true);
});

btnResetMeasure.addEventListener("click", () => {
  resetMeasurement();
  updateLabels();
  renderTable();
});

btnCSV.addEventListener("click", () => {
  const rows = [["t(s)", "theta(rad)", "omega(rad/s)", "alpha(rad/s^2)"]];
  for (const p of series) rows.push([p.t.toFixed(6), p.theta.toFixed(10), p.omega.toFixed(10), p.alpha.toFixed(10)]);
  const csv = rows.map(r => r.join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "pendulo_datos.csv";
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
});

// Sliders: when changed, update and reset
function bindResetOnInput(el, cb) {
  el.addEventListener("input", () => { cb(); updateLabels(); draw(); });
  el.addEventListener("change", () => { resetState(true); });
}

bindResetOnInput(Ls, () => {});
bindResetOnInput(ms, () => {});
bindResetOnInput(th0s, () => { theta = deg2rad(+th0s.value); omega = 0; });
bindResetOnInput(gs, () => {});
bindResetOnInput(betas, () => {});

nonlinearToggle.addEventListener("change", () => resetState(true));
dampingToggle.addEventListener("change", () => resetState(true));

// Main loop
function tick() {
  if (running && !dragging) {
    stepRK4();
    updateLabels();
    draw();
    drawPlots();
    renderTable();
  }
  requestAnimationFrame(tick);
}

// init
syncAllCanvasSizes();
resetState(true);
updateLabels();
renderTable();
tick();