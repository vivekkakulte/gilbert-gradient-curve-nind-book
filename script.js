const API_TUBING_CATALOG = [
  { od: '1.900"', id: 1.61, label: '1.900" OD (1.610" ID)' },
  { od: '2-3/8"', id: 1.995, label: '2-3/8" OD (1.995" ID)' },
  { od: '2-7/8"', id: 2.441, label: '2-7/8" OD (2.441" ID)' },
  { od: '3-1/2"', id: 2.992, label: '3-1/2" OD (2.992" ID)' },
  { od: '4-1/2"', id: 3.958, label: '4-1/2" OD (3.958" ID)' },
];

const REF_GLRS = [100, 200, 400, 600, 800, 1000, 1500, 2000, 3000];

let activeChartTab = "gilbert";
let userEditedGuess = false;
let lastSolvedState = null;

// Default Multi-Case Field Data Table Rows (each row has 1 blank target cell)
let batchRows = [
  {
    id: "Well Nind-01 (Benchmark)",
    depth: "8000",
    tubing: "1.995",
    rate: "600",
    glr: "400",
    pwh: "",
    pwf: "2200",
    solvedField: null,
    status: "Ready (Missing THP)",
  },
  {
    id: "Well Alpha-02",
    depth: "7500",
    tubing: "1.995",
    rate: "500",
    glr: "550",
    pwh: "220",
    pwf: "",
    solvedField: null,
    status: "Ready (Missing BHP)",
  },
  {
    id: "Well Beta-07 (Gas-Lift)",
    depth: "6500",
    tubing: "1.995",
    rate: "500",
    glr: "",
    pwh: "180",
    pwf: "1650",
    solvedField: null,
    status: "Ready (Missing GLR)",
  },
  {
    id: "Well Gamma-11",
    depth: "",
    tubing: "2.441",
    rate: "850",
    glr: "450",
    pwh: "240",
    pwf: "1950",
    solvedField: null,
    status: "Ready (Missing Depth)",
  },
];

function getGradient(p, q, glr, id) {
  const p_abs = Math.max(14.7, p + 14.7);
  const d_nom = id / 1.995;
  const q_ratio = Math.max(1.0, q) / 600.0;
  const glr_eff = Math.max(25.0, glr);
  const g_surface =
    0.12 * Math.pow(400.0 / glr_eff, 0.35) * Math.pow(1.0 / d_nom, 0.2);
  const g_deep = 0.266 * Math.pow(1.0 / d_nom, 0.1);
  const p_trans = 950.0 * Math.pow(glr_eff / 400.0, 0.5);
  const p_factor = Math.pow(Math.max(0, p) / (Math.max(0, p) + p_trans), 0.95);
  const grad_static = g_surface + (g_deep - g_surface) * p_factor;
  const v_rel =
    (q_ratio / Math.pow(d_nom, 2)) *
    (1.0 + (glr_eff * 14.7) / (5.615 * p_abs * 6.0));
  const grad_friction =
    0.0045 *
    Math.pow(1.0 / d_nom, 1.2) *
    Math.pow(v_rel, 1.4) *
    (p_abs / (p_abs + 400.0));
  return Math.min(0.45, Math.max(0.04, grad_static + grad_friction));
}

// Integrate H(P) from 0 to P using 4th-Order Runge-Kutta (step = 20 psi)
function computeEquivalentDepth(targetP, q, glr, id) {
  if (targetP <= 0) return 0;
  let p = 0;
  let h = 0;
  const step = 20;
  while (p < targetP) {
    const dp = Math.min(step, targetP - p);
    const k1 = 1.0 / getGradient(p, q, glr, id);
    const k2 = 1.0 / getGradient(p + 0.5 * dp, q, glr, id);
    const k3 = 1.0 / getGradient(p + 0.5 * dp, q, glr, id);
    const k4 = 1.0 / getGradient(p + dp, q, glr, id);
    h += (dp / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    p += dp;
  }
  return h;
}

// Find nearest standard API tubing OD for a continuous ID
function getNearestApiTubing(idInches) {
  let best = API_TUBING_CATALOG[0];
  let minDiff = Math.abs(idInches - best.id);
  for (const item of API_TUBING_CATALOG) {
    const diff = Math.abs(idInches - item.id);
    if (diff < minDiff) {
      minDiff = diff;
      best = item;
    }
  }
  return best;
}

function computeSmartInitialGuess(target, params) {
  const { depth, rate, glr, pwh, pwf } = params;
  if (target === "pwh") return Math.max(25, Math.round(pwf - 0.21 * depth));
  if (target === "pwf") return Math.round(pwh + 0.21 * depth);
  if (target === "depth") return Math.max(100, Math.round((pwf - pwh) / 0.21));
  if (target === "glr") return 500;
  if (target === "rate") return 500;
  if (target === "tubing") return 2.15;
  return 500;
}

function solveNewtonRaphson(
  target,
  params,
  tol = 1e-6,
  maxIter = 40,
  customX0 = null,
) {
  const { depth, id, rate, glr, pwh, pwf } = params;
  const history = [];

  const bounds = {
    pwh: { min: 0.0, max: 10000.0 },
    pwf: { min: 0.0, max: 15000.0 },
    depth: { min: 0.0, max: 30000.0 },
    glr: { min: 15.0, max: 15000.0 },
    rate: { min: 5.0, max: 25000.0 },
    tubing: { min: 0.75, max: 8.0 },
  };

  // Check non-flowing boundary when solving for Pwh
  if (target === "pwh") {
    const h2_check = computeEquivalentDepth(pwf, rate, glr, id);
    if (h2_check - depth <= 0) {
      const f0 = 0 - (h2_check - depth);
      const df0 = 1.0 / getGradient(0, rate, glr, id);
      history.push({ k: 1, x: 0, fx: f0, dfx: df0, dx: 0, relErr: 0 });
      return {
        root: 0,
        converged: false,
        nonFlowing: true,
        history,
        finalResidual: Math.abs(f0),
      };
    }
  }

  let x =
    customX0 !== null && !isNaN(customX0)
      ? customX0
      : computeSmartInitialGuess(target, params);
  x = Math.min(bounds[target].max, Math.max(bounds[target].min, x));

  // Root objective function f(x) = 0
  function evalF(val) {
    if (target === "pwh") {
      return (
        computeEquivalentDepth(val, rate, glr, id) -
        (computeEquivalentDepth(pwf, rate, glr, id) - depth)
      );
    }
    if (target === "pwf") {
      return (
        computeEquivalentDepth(val, rate, glr, id) -
        (computeEquivalentDepth(pwh, rate, glr, id) + depth)
      );
    }
    if (target === "depth") {
      return (
        computeEquivalentDepth(pwh, rate, glr, id) +
        val -
        computeEquivalentDepth(pwf, rate, glr, id)
      );
    }
    if (target === "glr") {
      return (
        computeEquivalentDepth(pwf, rate, val, id) -
        computeEquivalentDepth(pwh, rate, val, id) -
        depth
      );
    }
    if (target === "rate") {
      return (
        computeEquivalentDepth(pwf, val, glr, id) -
        computeEquivalentDepth(pwh, val, glr, id) -
        depth
      );
    }
    if (target === "tubing") {
      return (
        computeEquivalentDepth(pwf, rate, glr, val) -
        computeEquivalentDepth(pwh, rate, glr, val) -
        depth
      );
    }
    return 0;
  }

  // Exact Fundamental Theorem of Calculus derivative for Pwh/Pwf; damped central differences for GLR/Rate/ID
  function evalDf(val) {
    if (target === "pwh" || target === "pwf") {
      return 1.0 / getGradient(val, rate, glr, id);
    }
    if (target === "depth") {
      return 1.0;
    }
    const eps = Math.max(1e-4, 1e-4 * Math.abs(val));
    const xPlus = Math.min(bounds[target].max, val + eps);
    const xMinus = Math.max(bounds[target].min, val - eps);
    const span = xPlus - xMinus;
    if (span <= 0) return 1e-6;
    return (evalF(xPlus) - evalF(xMinus)) / span;
  }

  let converged = false;
  let finalResidual = 0;

  for (let iter = 1; iter <= maxIter; iter++) {
    const fx = evalF(x);
    const dfx = evalDf(x);
    const safeDfx = Math.abs(dfx) < 1e-12 ? (dfx >= 0 ? 1e-12 : -1e-12) : dfx;

    let dx = -fx / safeDfx;

    // Damped step guard for strongly nonlinear inversions (GLR, Rate, Tubing ID)
    if (target === "glr" || target === "rate" || target === "tubing") {
      const maxStep = 0.55 * Math.max(1.0, Math.abs(x));
      if (Math.abs(dx) > maxStep) {
        dx = Math.sign(dx) * maxStep;
      }
    }

    let xNext = x + dx;
    xNext = Math.min(bounds[target].max, Math.max(bounds[target].min, xNext));
    const actualDx = xNext - x;
    const relErr = Math.abs(actualDx) / Math.max(1e-8, Math.abs(xNext));
    const nextResidual = Math.abs(evalF(xNext));

    history.push({
      k: iter,
      x: xNext,
      fx: fx,
      dfx: dfx,
      dx: actualDx,
      relErr: relErr,
    });

    x = xNext;
    finalResidual = nextResidual;

    if (relErr <= tol && nextResidual <= Math.max(1e-3, tol * 100)) {
      converged = true;
      break;
    }
  }

  return {
    root: x,
    converged,
    nonFlowing: false,
    history,
    finalResidual,
  };
}

function updateMathFormulationCard(target) {
  const objEl = document.getElementById("mathObjectiveEq");
  const derEl = document.getElementById("mathDerivativeEq");

  if (target === "pwh") {
    objEl.innerHTML =
      "f(P<sub>wh</sub>) = H(P<sub>wh</sub>) &minus; [H(P<sub>wf</sub>) &minus; D] = 0";
    derEl.innerHTML =
      "f'(P<sub>wh</sub>) = dH/dP<sub>wh</sub> = 1 / (dp/dh)|<sub>Pwh</sub>";
  } else if (target === "pwf") {
    objEl.innerHTML =
      "f(P<sub>wf</sub>) = H(P<sub>wf</sub>) &minus; [H(P<sub>wh</sub>) + D] = 0";
    derEl.innerHTML =
      "f'(P<sub>wf</sub>) = dH/dP<sub>wf</sub> = 1 / (dp/dh)|<sub>Pwf</sub>";
  } else if (target === "depth") {
    objEl.innerHTML =
      "f(D) = [H(P<sub>wh</sub>) + D] &minus; H(P<sub>wf</sub>) = 0";
    derEl.innerHTML = "f'(D) = 1.0 (Exact Linear Depth Balance)";
  } else if (target === "glr") {
    objEl.innerHTML =
      "f(GLR) = [H(P<sub>wf</sub>, GLR) &minus; H(P<sub>wh</sub>, GLR)] &minus; D = 0";
    derEl.innerHTML =
      "f'(GLR) &approx; [f(GLR + &delta;) &minus; f(GLR &minus; &delta;)] / 2&delta;";
  } else if (target === "rate") {
    objEl.innerHTML =
      "f(q) = [H(P<sub>wf</sub>, q) &minus; H(P<sub>wh</sub>, q)] &minus; D = 0";
    derEl.innerHTML =
      "f'(q) &approx; [f(q + &delta;) &minus; f(q &minus; &delta;)] / 2&delta;";
  } else if (target === "tubing") {
    objEl.innerHTML =
      "f(d) = [H(P<sub>wf</sub>, d) &minus; H(P<sub>wh</sub>, d)] &minus; D = 0";
    derEl.innerHTML =
      "f'(d) &approx; [f(d + &delta;) &minus; f(d &minus; &delta;)] / 2&delta;";
  }
}

function handleTargetChange() {
  const target = document.getElementById("targetVar").value;
  const fields = ["depth", "tubing", "rate", "glr", "pwf", "pwh"];

  fields.forEach((f) => {
    const badge = document.getElementById(`badge-${f}`);
    if (f === "tubing") {
      const idInput = document.getElementById("tubingID");
      const sel = document.getElementById("tubingPreset");
      if (target === "tubing") {
        badge.classList.remove("hidden");
        idInput.setAttribute("readonly", "true");
        sel.disabled = true;
        idInput.className =
          "w-full bg-blue-50 border border-blue-400 rounded-lg px-3 py-2 pr-11 text-sm font-mono font-bold text-blue-900 focus:outline-none";
        sel.className =
          "w-full bg-blue-50 border border-blue-400 rounded-lg px-2.5 py-2 text-xs font-bold text-blue-900 opacity-85 cursor-not-allowed";
      } else {
        badge.classList.add("hidden");
        idInput.removeAttribute("readonly");
        sel.disabled = false;
        idInput.className =
          "w-full bg-white border border-slate-300 rounded-lg px-3 py-2 pr-11 text-sm font-mono text-slate-900 focus:outline-none focus:ring-2 focus:ring-blue-600";
        sel.className =
          "w-full bg-white border border-slate-300 rounded-lg px-2.5 py-2 text-xs font-medium text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-600";
      }
    } else {
      const input = document.getElementById(f);
      if (f === target) {
        badge.classList.remove("hidden");
        input.setAttribute("readonly", "true");
        input.className =
          "w-full bg-blue-50 border border-blue-400 rounded-lg px-3 py-2 pr-14 text-sm font-mono font-bold text-blue-900 focus:outline-none";
      } else {
        badge.classList.add("hidden");
        input.removeAttribute("readonly");
        input.className =
          "w-full bg-white border border-slate-300 rounded-lg px-3 py-2 pr-14 text-sm font-mono text-slate-900 focus:outline-none focus:ring-2 focus:ring-blue-600";
      }
    }
  });

  userEditedGuess = false;
  updateMathFormulationCard(target);
  calculate(false);
}

function syncTubingFromDropdown() {
  const val = document.getElementById("tubingPreset").value;
  if (val !== "custom") {
    document.getElementById("tubingID").value = parseFloat(val).toFixed(3);
  }
  updateNearestApiDisplay(
    parseFloat(document.getElementById("tubingID").value) || 1.995,
  );
  onInputModified();
}

function syncTubingFromInput() {
  const idVal = parseFloat(document.getElementById("tubingID").value) || 1.995;
  const match = API_TUBING_CATALOG.find(
    (item) => Math.abs(item.id - idVal) < 0.002,
  );
  document.getElementById("tubingPreset").value = match
    ? match.id.toFixed(3)
    : "custom";
  updateNearestApiDisplay(idVal);
  onInputModified();
}

function updateNearestApiDisplay(idVal) {
  const nearest = getNearestApiTubing(idVal);
  document.getElementById("nearestApiText").textContent = nearest.label;
}

function resetAutoGuess() {
  userEditedGuess = false;
  calculate(true);
  showToast("Initial guess (x₀) reset to automatic estimate");
}

function onInputModified() {
  calculate(false);
}

function calculate(notify = false) {
  const target = document.getElementById("targetVar").value;
  const tolRaw = parseFloat(document.getElementById("nrTolerance").value);
  const tol = !isNaN(tolRaw) && tolRaw > 0 ? tolRaw : 1e-6;

  let params = {
    depth: Math.max(
      0,
      parseFloat(document.getElementById("depth").value) || 8000,
    ),
    id: Math.max(
      0.5,
      parseFloat(document.getElementById("tubingID").value) || 1.995,
    ),
    rate: Math.max(1, parseFloat(document.getElementById("rate").value) || 600),
    glr: Math.max(15, parseFloat(document.getElementById("glr").value) || 400),
    pwf: Math.max(0, parseFloat(document.getElementById("pwf").value) || 2200),
    pwh: Math.max(0, parseFloat(document.getElementById("pwh").value) || 530),
  };

  const guessInput = document.getElementById("nrGuess");
  if (!userEditedGuess) {
    const autoX0 = computeSmartInitialGuess(target, params);
    guessInput.value = autoX0;
  }
  const x0 = parseFloat(guessInput.value);

  const nrResult = solveNewtonRaphson(target, params, tol, 40, x0);
  const solved = nrResult.root;

  // Write solved value to the target input field
  if (target === "pwh") {
    params.pwh = solved;
    document.getElementById("pwh").value = solved.toFixed(1);
  } else if (target === "pwf") {
    params.pwf = solved;
    document.getElementById("pwf").value = solved.toFixed(1);
  } else if (target === "depth") {
    params.depth = solved;
    document.getElementById("depth").value = Math.round(solved);
  } else if (target === "glr") {
    params.glr = solved;
    document.getElementById("glr").value = Math.round(solved);
  } else if (target === "rate") {
    params.rate = solved;
    document.getElementById("rate").value = Math.round(solved);
  } else if (target === "tubing") {
    params.id = solved;
    document.getElementById("tubingID").value = solved.toFixed(3);
    const nearest = getNearestApiTubing(solved);
    document.getElementById("tubingPreset").value = nearest.id.toFixed(3);
    updateNearestApiDisplay(solved);
  }

  const h1 = computeEquivalentDepth(
    params.pwh,
    params.rate,
    params.glr,
    params.id,
  );
  const h2 = h1 + params.depth;
  const deltaP = params.pwf - params.pwh;
  const avgGrad = params.depth > 0 ? deltaP / params.depth : 0;

  lastSolvedState = {
    target,
    ...params,
    h1,
    h2,
    deltaP,
    avgGrad,
    tol,
    nrResult,
  };

  // Update Summary Metrics
  document.getElementById("resDeltaP").textContent = `${deltaP.toFixed(1)} psi`;
  document.getElementById("resAvgGrad").textContent =
    `${avgGrad.toFixed(4)} psi/ft`;
  document.getElementById("resH1").textContent =
    `${Math.round(h1).toLocaleString()} ft`;
  document.getElementById("resH2").textContent =
    `${Math.round(h2).toLocaleString()} ft`;
  document.getElementById("legendGLR").textContent = Math.round(params.glr);
  document.getElementById("graphicalBalanceText").textContent =
    `H₂ (${Math.round(h2).toLocaleString()} ft) − H₁ (${Math.round(h1).toLocaleString()} ft) = ${Math.round(params.depth).toLocaleString()} ft`;

  renderConvergenceHistory(nrResult, target);
  renderTraverseTable();
  runVarianceAnalysis();
  drawCanvas();

  if (notify) {
    showToast(
      `Newton-Raphson converged in ${nrResult.history.length} iteration(s)`,
    );
  }
}

function renderConvergenceHistory(nrResult, target) {
  const tbody = document.getElementById("nrHistoryBody");
  const badge = document.getElementById("nrConvBadge");
  const targetLabelEl = document.getElementById("nrHeaderTarget");
  const resLabelEl = document.getElementById("nrHeaderResidual");

  const targetNames = {
    pwh: "THP / Pwh (psig)",
    pwf: "BHP / Pwf (psig)",
    depth: "Wellbore Depth D (ft)",
    glr: "Gas-Liquid Ratio GLR (scf/bbl)",
    rate: "Liquid Rate q (bbl/day)",
    tubing: "Tubing ID (in)",
  };
  targetLabelEl.textContent = targetNames[target] || target;
  resLabelEl.textContent = nrResult.finalResidual.toExponential(2);

  if (nrResult.nonFlowing) {
    badge.innerHTML =
      '<i class="fa-solid fa-triangle-exclamation mr-1"></i> Non-Flowing (THP = 0)';
    badge.className =
      "px-3 py-1 rounded-full bg-red-50 text-red-700 border border-red-200 font-bold";
  } else if (nrResult.converged) {
    badge.innerHTML = `<i class="fa-solid fa-circle-check mr-1"></i> Converged in ${nrResult.history.length} iterations`;
    badge.className =
      "px-3 py-1 rounded-full bg-emerald-50 text-emerald-700 border border-emerald-200 font-bold";
  } else {
    badge.innerHTML = `<i class="fa-solid fa-clock mr-1"></i> Stopped at ${nrResult.history.length} iterations`;
    badge.className =
      "px-3 py-1 rounded-full bg-amber-50 text-amber-700 border border-amber-200 font-bold";
  }

  tbody.innerHTML = "";
  nrResult.history.forEach((row, idx) => {
    const isLast = idx === nrResult.history.length - 1;
    const tr = document.createElement("tr");
    tr.className = isLast
      ? "bg-emerald-50/70 font-semibold text-slate-900"
      : idx % 2 === 0
        ? "bg-white"
        : "bg-slate-50/60";

    tr.innerHTML = `
          <td class="py-2 px-4 font-bold text-blue-900">${row.k}</td>
          <td class="py-2 px-4 text-blue-700 font-semibold">${row.x.toFixed(4)}</td>
          <td class="py-2 px-4">${row.fx.toExponential(3)}</td>
          <td class="py-2 px-4 text-slate-600">${row.dfx.toFixed(4)}</td>
          <td class="py-2 px-4">${row.dx >= 0 ? "+" : ""}${row.dx.toExponential(3)}</td>
          <td class="py-2 px-4 ${isLast ? "text-emerald-700 font-bold" : "text-slate-600"}">${(row.relErr * 100).toFixed(5)}%</td>
        `;
    tbody.appendChild(tr);
  });
}

// Invert H(P) = targetH via fast Newton-Raphson for depth stations
function invertEquivalentDepth(targetH, q, glr, id) {
  if (targetH <= 0) return 0;
  let p = Math.min(9000, targetH * 0.21);
  for (let i = 0; i < 15; i++) {
    const f = computeEquivalentDepth(p, q, glr, id) - targetH;
    if (Math.abs(f) < 1e-4) break;
    const df = 1.0 / getGradient(p, q, glr, id);
    p = Math.max(0, p - f / df);
  }
  return p;
}

function switchChartTab(tab) {
  activeChartTab = tab;
  const btnG = document.getElementById("btnTabGilbert");
  const btnT = document.getElementById("btnTabTraverse");
  const refLegend = document.getElementById("legendRefWrap");

  if (tab === "gilbert") {
    btnG.className =
      "px-3 py-1.5 text-xs font-semibold rounded-md bg-blue-600 text-white shadow-sm transition";
    btnT.className =
      "px-3 py-1.5 text-xs font-semibold rounded-md text-slate-600 hover:text-slate-900 transition";
    refLegend.classList.remove("hidden");
  } else {
    btnT.className =
      "px-3 py-1.5 text-xs font-semibold rounded-md bg-blue-600 text-white shadow-sm transition";
    btnG.className =
      "px-3 py-1.5 text-xs font-semibold rounded-md text-slate-600 hover:text-slate-900 transition";
    refLegend.classList.add("hidden");
  }
  drawCanvas();
}

function drawCanvas() {
  const canvas = document.getElementById("mainCanvas");
  if (!canvas || !lastSolvedState) return;
  const ctx = canvas.getContext("2d");
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;

  canvas.width = rect.width * dpr;
  canvas.height = rect.height * dpr;
  ctx.scale(dpr, dpr);

  const w = rect.width;
  const h = rect.height;
  const pad = { top: 36, right: 36, bottom: 48, left: 68 };
  const plotW = w - pad.left - pad.right;
  const plotH = h - pad.top - pad.bottom;

  // Pure white engineering paper background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, w, h);

  if (activeChartTab === "gilbert") {
    drawGilbertGraphPaper(ctx, pad, plotW, plotH, w, h);
  } else {
    drawTraverseGraphPaper(ctx, pad, plotW, plotH, w, h);
  }
}

function drawGilbertGraphPaper(ctx, pad, pw, ph, totalW, totalH) {
  const s = lastSolvedState;
  const maxP = Math.max(3500, Math.ceil((s.pwf + 400) / 500) * 500);
  const maxH = Math.max(14000, Math.ceil((s.h2 + 1500) / 2000) * 2000);

  const toX = (p) => pad.left + (p / maxP) * pw;
  const toY = (d) => pad.top + (d / maxH) * ph; // Inverted depth axis: H=0 at top

  // Crisp #e2e8f0 graph paper gridlines
  ctx.strokeStyle = "#e2e8f0";
  ctx.lineWidth = 1;
  ctx.fillStyle = "#475569";
  ctx.font = "10px JetBrains Mono, monospace";

  const pStep = maxP > 4500 ? 1000 : 500;
  ctx.textAlign = "center";
  for (let p = 0; p <= maxP; p += pStep) {
    const x = toX(p);
    ctx.beginPath();
    ctx.moveTo(x, pad.top);
    ctx.lineTo(x, pad.top + ph);
    ctx.stroke();
    ctx.fillText(p, x, pad.top + ph + 17);
  }

  const hStep = 2000;
  ctx.textAlign = "right";
  for (let d = 0; d <= maxH; d += hStep) {
    const y = toY(d);
    ctx.beginPath();
    ctx.moveTo(pad.left, y);
    ctx.lineTo(pad.left + pw, y);
    ctx.stroke();
    ctx.fillText(d.toLocaleString(), pad.left - 10, y + 3);
  }

  // Plot frame border
  ctx.strokeStyle = "#cbd5e1";
  ctx.lineWidth = 1.2;
  ctx.strokeRect(pad.left, pad.top, pw, ph);

  // Family of reference GLR curves in soft gray (#94a3b8) with inline curve labels
  REF_GLRS.forEach((refGlr) => {
    if (Math.abs(refGlr - s.glr) < 5) return;
    ctx.beginPath();
    ctx.strokeStyle = "#94a3b8";
    ctx.lineWidth = 1.2;
    let started = false;
    let labelPt = null;

    for (let p = 0; p <= maxP; p += 40) {
      const eqH = computeEquivalentDepth(p, s.rate, refGlr, s.id);
      if (eqH > maxH) break;
      const x = toX(p);
      const y = toY(eqH);
      if (!started) {
        ctx.moveTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
      if (eqH <= maxH * 0.94 && x <= pad.left + pw - 38) {
        labelPt = { x, y };
      }
    }
    ctx.stroke();

    if (labelPt) {
      ctx.fillStyle = "#64748b";
      ctx.font = "600 9px JetBrains Mono, monospace";
      ctx.textAlign = "left";
      ctx.fillText(`${refGlr}`, labelPt.x + 4, labelPt.y - 2);
    }
  });

  // Active GLR curve in bold Royal Blue (#2563eb)
  ctx.beginPath();
  ctx.strokeStyle = "#2563eb";
  ctx.lineWidth = 3;
  let activeStarted = false;
  for (let p = 0; p <= maxP; p += 20) {
    const eqH = computeEquivalentDepth(p, s.rate, s.glr, s.id);
    if (eqH > maxH) break;
    const x = toX(p);
    const y = toY(eqH);
    if (!activeStarted) {
      ctx.moveTo(x, y);
      activeStarted = true;
    } else {
      ctx.lineTo(x, y);
    }
  }
  ctx.stroke();

  // Coordinates for Point 1 (Wellhead: Pwh, H1) and Point 2 (Sandface: Pwf, H2)
  const x1 = toX(s.pwh);
  const y1 = toY(s.h1);
  const x2 = toX(s.pwf);
  const y2 = toY(s.h2);

  // Dashed orthogonal projections to both axes
  ctx.setLineDash([4, 4]);
  ctx.lineWidth = 1.4;

  // Crimson projections for Wellhead (Pwh, H1)
  ctx.strokeStyle = "#dc2626";
  ctx.beginPath();
  ctx.moveTo(pad.left, y1);
  ctx.lineTo(x1, y1);
  ctx.lineTo(x1, pad.top + ph);
  ctx.stroke();

  // Emerald projections for Sandface (Pwf, H2)
  ctx.strokeStyle = "#059669";
  ctx.beginPath();
  ctx.moveTo(pad.left, y2);
  ctx.lineTo(x2, y2);
  ctx.lineTo(x2, pad.top + ph);
  ctx.stroke();
  ctx.setLineDash([]);

  // Vertical dimension bracket labeled "ΔH = D"
  const bracketX = pad.left + pw - 24;
  ctx.strokeStyle = "#1e3a8a";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(bracketX - 6, y1);
  ctx.lineTo(bracketX, y1);
  ctx.lineTo(bracketX, y2);
  ctx.lineTo(bracketX - 6, y2);
  ctx.stroke();

  ctx.fillStyle = "#1e3a8a";
  ctx.font = "bold 11px JetBrains Mono, monospace";
  ctx.textAlign = "right";
  ctx.fillText(
    `ΔH = D (${Math.round(s.depth).toLocaleString()} ft)`,
    bracketX - 10,
    (y1 + y2) / 2 + 4,
  );

  // Markers
  drawPaperMarker(ctx, x1, y1, "#dc2626");
  drawPaperMarker(ctx, x2, y2, "#059669");

  // Clean white callout labels
  drawPaperTag(
    ctx,
    x1 + 10,
    y1 - 8,
    `Pwh = ${s.pwh.toFixed(1)} psig (H₁ = ${Math.round(s.h1)} ft)`,
    "#dc2626",
    "#fef2f2",
  );
  drawPaperTag(
    ctx,
    Math.min(x2 + 10, bracketX - 210),
    y2 + 16,
    `Pwf = ${s.pwf.toFixed(1)} psig (H₂ = ${Math.round(s.h2)} ft)`,
    "#059669",
    "#ecfdf5",
  );

  // Axis Titles
  ctx.fillStyle = "#1e3a8a";
  ctx.font = "bold 11px Inter, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(
    "Flowing Tubing Pressure P (psig)",
    pad.left + pw / 2,
    totalH - 12,
  );

  ctx.save();
  ctx.translate(18, pad.top + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Gilbert Equivalent Depth H (ft) [H = 0 at Top]", 0, 0);
  ctx.restore();
}

function drawTraverseGraphPaper(ctx, pad, pw, ph, totalW, totalH) {
  const s = lastSolvedState;
  const staticTD = s.pwh + s.depth * 0.368;
  const maxP = Math.max(
    3000,
    Math.ceil((Math.max(s.pwf, staticTD) + 400) / 500) * 500,
  );
  const maxD = Math.max(1000, s.depth);

  const toX = (p) => pad.left + (p / maxP) * pw;
  const toY = (z) => pad.top + (z / maxD) * ph;

  // Gridlines
  ctx.strokeStyle = "#e2e8f0";
  ctx.lineWidth = 1;
  ctx.fillStyle = "#475569";
  ctx.font = "10px JetBrains Mono, monospace";

  const pStep = maxP > 4500 ? 1000 : 500;
  ctx.textAlign = "center";
  for (let p = 0; p <= maxP; p += pStep) {
    const x = toX(p);
    ctx.beginPath();
    ctx.moveTo(x, pad.top);
    ctx.lineTo(x, pad.top + ph);
    ctx.stroke();
    ctx.fillText(p, x, pad.top + ph + 17);
  }

  ctx.textAlign = "right";
  for (let i = 0; i <= 5; i++) {
    const z = (i / 5) * maxD;
    const y = toY(z);
    ctx.beginPath();
    ctx.moveTo(pad.left, y);
    ctx.lineTo(pad.left + pw, y);
    ctx.stroke();
    ctx.fillText(Math.round(z).toLocaleString(), pad.left - 10, y + 3);
  }

  ctx.strokeStyle = "#cbd5e1";
  ctx.lineWidth = 1.2;
  ctx.strokeRect(pad.left, pad.top, pw, ph);

  // Static liquid gradient reference (0.368 psi/ft)
  ctx.setLineDash([5, 5]);
  ctx.strokeStyle = "#94a3b8";
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ctx.moveTo(toX(s.pwh), toY(0));
  ctx.lineTo(toX(staticTD), toY(s.depth));
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "#64748b";
  ctx.font = "600 10px JetBrains Mono, monospace";
  ctx.textAlign = "right";
  ctx.fillText(
    `Static Liquid Gradient (0.368 psi/ft = ${Math.round(staticTD)} psig)`,
    toX(staticTD) - 8,
    toY(s.depth) - 10,
  );

  // Multiphase Flowing Pressure Traverse in Royal Blue (#2563eb)
  ctx.beginPath();
  ctx.strokeStyle = "#2563eb";
  ctx.lineWidth = 3;
  const steps = 40;
  for (let i = 0; i <= steps; i++) {
    const z = (i / steps) * s.depth;
    const p =
      i === 0
        ? s.pwh
        : i === steps
          ? s.pwf
          : invertEquivalentDepth(s.h1 + z, s.rate, s.glr, s.id);
    const x = toX(p);
    const y = toY(z);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  const x0 = toX(s.pwh);
  const y0 = toY(0);
  const xTD = toX(s.pwf);
  const yTD = toY(s.depth);

  drawPaperMarker(ctx, x0, y0, "#dc2626");
  drawPaperMarker(ctx, xTD, yTD, "#059669");

  drawPaperTag(
    ctx,
    x0 + 10,
    y0 + 16,
    `Wellhead Pwh = ${s.pwh.toFixed(1)} psig @ 0 ft`,
    "#dc2626",
    "#fef2f2",
  );
  drawPaperTag(
    ctx,
    Math.max(pad.left + 10, xTD - 205),
    yTD - 12,
    `Sandface Pwf = ${s.pwf.toFixed(1)} psig @ ${Math.round(s.depth)} ft`,
    "#059669",
    "#ecfdf5",
  );

  ctx.fillStyle = "#1e3a8a";
  ctx.font = "bold 11px Inter, sans-serif";
  ctx.textAlign = "center";
  ctx.fillText(
    "Flowing Wellbore Pressure P(z) (psig)",
    pad.left + pw / 2,
    totalH - 12,
  );

  ctx.save();
  ctx.translate(18, pad.top + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("True Vertical Well Depth z (ft)", 0, 0);
  ctx.restore();
}

function drawPaperMarker(ctx, x, y, color) {
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(x, y, 5.5, 0, Math.PI * 2);
  ctx.fill();
  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 2;
  ctx.stroke();
}

function drawPaperTag(ctx, x, y, text, borderColor, bgColor) {
  ctx.font = "bold 10px JetBrains Mono, monospace";
  const w = ctx.measureText(text).width;
  ctx.fillStyle = bgColor;
  ctx.strokeStyle = borderColor;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.roundRect(x - 5, y - 10, w + 10, 16, 4);
  ctx.fill();
  ctx.stroke();
  ctx.fillStyle = borderColor;
  ctx.textAlign = "left";
  ctx.fillText(text, x, y + 1);
}

// Hover & Touch Coordinate Tooltip on Graph Paper Canvas
const canvasEl = document.getElementById("mainCanvas");
const tooltipEl = document.getElementById("canvasTooltip");

function handleCanvasPointer(clientX, clientY) {
  if (!lastSolvedState) return;
  const rect = canvasEl.getBoundingClientRect();
  const mx = clientX - rect.left;
  const my = clientY - rect.top;
  const pad = { top: 36, right: 36, bottom: 48, left: 68 };
  const pw = rect.width - pad.left - pad.right;
  const ph = rect.height - pad.top - pad.bottom;

  if (
    mx >= pad.left &&
    mx <= pad.left + pw &&
    my >= pad.top &&
    my <= pad.top + ph
  ) {
    const s = lastSolvedState;
    tooltipEl.classList.remove("hidden");
    const leftPos = mx + 185 > rect.width ? mx - 195 : mx + 14;
    const topPos = my + 80 > rect.height ? my - 75 : my + 14;
    tooltipEl.style.left = `${leftPos}px`;
    tooltipEl.style.top = `${topPos}px`;

    if (activeChartTab === "gilbert") {
      const maxP = Math.max(3500, Math.ceil((s.pwf + 400) / 500) * 500);
      const maxH = Math.max(14000, Math.ceil((s.h2 + 1500) / 2000) * 2000);
      const p = ((mx - pad.left) / pw) * maxP;
      const h = ((my - pad.top) / ph) * maxH;
      const grad = getGradient(p, s.rate, s.glr, s.id);
      tooltipEl.innerHTML = `
            <div class="text-blue-800 font-bold mb-0.5">Gilbert Graph Coordinates</div>
            <div>Pressure P: <span class="font-bold text-slate-900">${Math.round(p)} psig</span></div>
            <div>Eq. Depth H: <span class="font-bold text-slate-900">${Math.round(h).toLocaleString()} ft</span></div>
            <div>Local dp/dh: <span class="font-bold text-blue-600">${grad.toFixed(4)} psi/ft</span></div>
          `;
    } else {
      const staticTD = s.pwh + s.depth * 0.368;
      const maxP = Math.max(
        3000,
        Math.ceil((Math.max(s.pwf, staticTD) + 400) / 500) * 500,
      );
      const p = ((mx - pad.left) / pw) * maxP;
      const z = ((my - pad.top) / ph) * Math.max(1000, s.depth);
      const grad = getGradient(p, s.rate, s.glr, s.id);
      tooltipEl.innerHTML = `
            <div class="text-blue-800 font-bold mb-0.5">Wellbore Traverse Cursor</div>
            <div>Depth z: <span class="font-bold text-slate-900">${Math.round(z).toLocaleString()} ft</span></div>
            <div>Pressure P: <span class="font-bold text-slate-900">${Math.round(p)} psig</span></div>
            <div>Local dp/dh: <span class="font-bold text-blue-600">${grad.toFixed(4)} psi/ft</span></div>
          `;
    }
  } else {
    tooltipEl.classList.add("hidden");
  }
}

canvasEl.addEventListener("mousemove", (e) =>
  handleCanvasPointer(e.clientX, e.clientY),
);
canvasEl.addEventListener("mouseleave", () =>
  tooltipEl.classList.add("hidden"),
);
canvasEl.addEventListener(
  "touchstart",
  (e) => {
    if (e.touches.length > 0)
      handleCanvasPointer(e.touches[0].clientX, e.touches[0].clientY);
  },
  { passive: true },
);
canvasEl.addEventListener(
  "touchmove",
  (e) => {
    if (e.touches.length > 0)
      handleCanvasPointer(e.touches[0].clientX, e.touches[0].clientY);
  },
  { passive: true },
);
canvasEl.addEventListener("touchend", () => tooltipEl.classList.add("hidden"));

function switchBottomTab(tab) {
  const tabs = ["variance", "batch", "traverse"];
  tabs.forEach((t) => {
    const btn = document.getElementById(
      `subTab${t.charAt(0).toUpperCase() + t.slice(1)}`,
    );
    const panel = document.getElementById(
      `panel${t.charAt(0).toUpperCase() + t.slice(1)}`,
    );
    if (t === tab) {
      btn.className =
        "px-3 py-1.5 text-xs font-semibold rounded-md bg-blue-600 text-white shadow-sm transition";
      panel.classList.remove("hidden");
    } else {
      btn.className =
        "px-3 py-1.5 text-xs font-semibold rounded-md text-slate-600 hover:text-slate-900 transition";
      panel.classList.add("hidden");
    }
  });
}

// Parameter Variance (±%) Sensitivity Analysis
function runVarianceAnalysis() {
  if (!lastSolvedState) return;
  const s = lastSolvedState;
  const target = s.target;
  let varParam = document.getElementById("varInputParam").value;
  const pct = parseFloat(document.getElementById("varPercentage").value) || 10;
  const autoNote = document.getElementById("varAutoSwitchNote");

  if (varParam === target) {
    varParam = target === "glr" ? "rate" : "glr";
    document.getElementById("varInputParam").value = varParam;
    autoNote.classList.remove("hidden");
  } else {
    autoNote.classList.add("hidden");
  }

  const paramKey = varParam === "tubing" ? "id" : varParam;
  const baseParams = {
    depth: s.depth,
    id: s.id,
    rate: s.rate,
    glr: s.glr,
    pwf: s.pwf,
    pwh: s.pwh,
  };

  const scenarios = [
    {
      label: `Low (-${pct}%)`,
      factor: 1 - pct / 100,
      badgeClass: "text-red-700 bg-red-50 border-red-200",
    },
    {
      label: "Nominal (0%)",
      factor: 1.0,
      badgeClass: "text-blue-800 bg-blue-50 border-blue-200",
    },
    {
      label: `High (+${pct}%)`,
      factor: 1 + pct / 100,
      badgeClass: "text-emerald-700 bg-emerald-50 border-emerald-200",
    },
  ];

  const tbody = document.getElementById("varianceTableBody");
  tbody.innerHTML = "";
  const solvedOutputs = [];

  scenarios.forEach((sc, idx) => {
    const runParams = {
      ...baseParams,
      [paramKey]: baseParams[paramKey] * sc.factor,
    };
    const res =
      idx === 1
        ? { root: s[target === "tubing" ? "id" : target] }
        : solveNewtonRaphson(target, runParams, s.tol, 35, null);

    const outVal = res.root;
    solvedOutputs.push(outVal);

    const caseFull = {
      ...runParams,
      [target === "tubing" ? "id" : target]: outVal,
    };
    const h1 = computeEquivalentDepth(
      caseFull.pwh,
      caseFull.rate,
      caseFull.glr,
      caseFull.id,
    );
    const h2 = h1 + caseFull.depth;
    const dP = caseFull.pwf - caseFull.pwh;

    const fmtParam = (k, v) => {
      if (k === "glr") return `${Math.round(v)} scf/bbl`;
      if (k === "rate") return `${Math.round(v)} bbl/day`;
      if (k === "depth") return `${Math.round(v).toLocaleString()} ft`;
      if (k === "pwf" || k === "pwh") return `${v.toFixed(1)} psig`;
      return `${v.toFixed(3)} in`;
    };

    const tr = document.createElement("tr");
    tr.className = idx === 1 ? "bg-blue-50/40 font-semibold" : "bg-white";
    tr.innerHTML = `
          <td class="py-2.5 px-4">
            <span class="px-2 py-0.5 rounded border text-[11px] font-bold ${sc.badgeClass}">${sc.label}</span>
          </td>
          <td class="py-2.5 px-4 text-slate-700">${fmtParam(varParam, runParams[paramKey])}</td>
          <td class="py-2.5 px-4 font-bold text-blue-900">${fmtParam(target, outVal)}</td>
          <td class="py-2.5 px-4 text-red-600">${Math.round(h1).toLocaleString()} ft</td>
          <td class="py-2.5 px-4 text-emerald-700">${Math.round(h2).toLocaleString()} ft</td>
          <td class="py-2.5 px-4 text-slate-800">${dP.toFixed(1)} psi</td>
        `;
    tbody.appendChild(tr);
  });

  const minOut = Math.min(...solvedOutputs);
  const maxOut = Math.max(...solvedOutputs);
  const diff = maxOut - minOut;
  const unitStr =
    target === "tubing"
      ? "in"
      : target === "pwh" || target === "pwf"
        ? "psi"
        : target === "depth"
          ? "ft"
          : target === "glr"
            ? "scf/bbl"
            : "bbl/d";
  const fmtVal =
    target === "tubing" ? (v) => v.toFixed(3) : (v) => v.toFixed(1);
  document.getElementById("varTotalSpread").textContent =
    `${fmtVal(diff)} ${unitStr} (${fmtVal(minOut)} to ${fmtVal(maxOut)} ${unitStr})`;
}

// Discrete Depth Traverse Profile Table
function renderTraverseTable() {
  const s = lastSolvedState;
  const tbody = document.getElementById("traverseTableBody");
  tbody.innerHTML = "";
  const steps = 10;
  const dz = s.depth / steps;

  for (let i = 0; i <= steps; i++) {
    const z = i * dz;
    const eqH = s.h1 + z;
    const p =
      i === 0
        ? s.pwh
        : i === steps
          ? s.pwf
          : invertEquivalentDepth(eqH, s.rate, s.glr, s.id);
    const grad = getGradient(p, s.rate, s.glr, s.id);

    const tr = document.createElement("tr");
    tr.className = i % 2 === 0 ? "bg-white" : "bg-slate-50/70";
    const stationLabel =
      i === 0 ? "Wellhead (0)" : i === steps ? "Sandface (TD)" : `Station ${i}`;
    const pressColor =
      i === 0
        ? "text-red-600 font-bold"
        : i === steps
          ? "text-emerald-700 font-bold"
          : "text-blue-800 font-semibold";

    tr.innerHTML = `
          <td class="py-2 px-4 text-slate-500">${stationLabel}</td>
          <td class="py-2 px-4 font-semibold text-slate-900">${Math.round(z).toLocaleString()}</td>
          <td class="py-2 px-4 ${pressColor}">${p.toFixed(1)}</td>
          <td class="py-2 px-4 text-slate-700">${grad.toFixed(4)}</td>
          <td class="py-2 px-4 text-slate-600">${Math.round(eqH).toLocaleString()}</td>
        `;
    tbody.appendChild(tr);
  }
}

// Multi-Case Field Data Table (Batch Solver)
function renderBatchTable() {
  const tbody = document.getElementById("batchTableBody");
  tbody.innerHTML = "";

  batchRows.forEach((row, idx) => {
    const tr = document.createElement("tr");
    tr.className = idx % 2 === 0 ? "bg-white" : "bg-slate-50/70";

    const inputCls = (field) =>
      row.solvedField === field
        ? "w-24 bg-blue-50 border border-blue-400 text-blue-900 font-bold rounded-md px-2 py-1 text-xs focus:outline-none"
        : "w-24 bg-white border border-slate-300 text-slate-900 rounded-md px-2 py-1 text-xs focus:outline-none focus:ring-1 focus:ring-blue-600";

    tr.innerHTML = `
          <td class="py-2 px-3">
            <input type="text" value="${row.id}" onchange="updateBatchCell(${idx}, 'id', this.value)" class="w-40 bg-white border border-slate-300 text-slate-800 rounded-md px-2 py-1 text-xs focus:outline-none focus:ring-1 focus:ring-blue-600">
          </td>
          <td class="py-2 px-3">
            <input type="number" placeholder="Blank" value="${row.depth}" onchange="updateBatchCell(${idx}, 'depth', this.value)" class="${inputCls("depth")}">
          </td>
          <td class="py-2 px-3">
            <input type="number" step="0.005" placeholder="Blank" value="${row.tubing}" onchange="updateBatchCell(${idx}, 'tubing', this.value)" class="${inputCls("tubing")}">
            ${row.apiNote ? `<div class="text-[10px] text-blue-700 font-semibold mt-0.5">${row.apiNote}</div>` : ""}
          </td>
          <td class="py-2 px-3">
            <input type="number" placeholder="Blank" value="${row.rate}" onchange="updateBatchCell(${idx}, 'rate', this.value)" class="${inputCls("rate")}">
          </td>
          <td class="py-2 px-3">
            <input type="number" placeholder="Blank" value="${row.glr}" onchange="updateBatchCell(${idx}, 'glr', this.value)" class="${inputCls("glr")}">
          </td>
          <td class="py-2 px-3">
            <input type="number" placeholder="Blank" value="${row.pwh}" onchange="updateBatchCell(${idx}, 'pwh', this.value)" class="${inputCls("pwh")}">
          </td>
          <td class="py-2 px-3">
            <input type="number" placeholder="Blank" value="${row.pwf}" onchange="updateBatchCell(${idx}, 'pwf', this.value)" class="${inputCls("pwf")}">
          </td>
          <td class="py-2 px-3">
            <span class="text-[11px] ${row.solvedField ? "text-emerald-700 font-bold" : "text-slate-500"}">${row.status}</span>
          </td>
          <td class="py-2 px-3 text-center whitespace-nowrap">
            <button onclick="loadBatchRow(${idx})" class="px-2 py-1 text-[11px] font-semibold bg-blue-50 hover:bg-blue-100 text-blue-700 border border-blue-200 rounded mr-1 transition">
              Load
            </button>
            <button onclick="deleteBatchRow(${idx})" class="px-2 py-1 text-[11px] font-semibold bg-red-50 hover:bg-red-100 text-red-600 border border-red-200 rounded transition">
              <i class="fa-solid fa-trash"></i>
            </button>
          </td>
        `;
    tbody.appendChild(tr);
  });
}

function updateBatchCell(idx, field, val) {
  batchRows[idx][field] = val.trim();
  if (batchRows[idx].solvedField === field) {
    batchRows[idx].solvedField = null;
    batchRows[idx].apiNote = null;
  }
  batchRows[idx].status = "Edited";
}

function addBatchRow() {
  batchRows.push({
    id: `Well Case #${batchRows.length + 1}`,
    depth: "7500",
    tubing: "1.995",
    rate: "500",
    glr: "450",
    pwh: "",
    pwf: "1950",
    solvedField: null,
    status: "Ready (Missing THP)",
  });
  renderBatchTable();
}

function deleteBatchRow(idx) {
  if (batchRows.length <= 1) return;
  batchRows.splice(idx, 1);
  renderBatchTable();
}

function clearBatchSolved() {
  batchRows.forEach((r) => {
    if (r.solvedField) {
      r[r.solvedField] = "";
      r.solvedField = null;
      r.apiNote = null;
      r.status = "Cleared";
    }
  });
  renderBatchTable();
  showToast("Cleared solved batch cells");
}

function solveAllBatchRows() {
  const tolRaw = parseFloat(document.getElementById("nrTolerance").value);
  const tol = !isNaN(tolRaw) && tolRaw > 0 ? tolRaw : 1e-6;
  const fields = ["depth", "tubing", "rate", "glr", "pwh", "pwf"];
  let solvedCount = 0;
  let skippedCount = 0;

  batchRows.forEach((row) => {
    if (row.solvedField) {
      row[row.solvedField] = "";
      row.solvedField = null;
      row.apiNote = null;
    }
    const blanks = fields.filter(
      (f) => row[f] === "" || isNaN(parseFloat(row[f])),
    );
    if (blanks.length === 0) {
      row.status = "Complete (6/6)";
      return;
    }
    if (blanks.length > 1) {
      row.status = `Leave only 1 blank (${blanks.length} blank)`;
      skippedCount++;
      return;
    }

    const target = blanks[0];
    const params = {
      depth: parseFloat(row.depth) || 7000,
      id: parseFloat(row.tubing) || 1.995,
      rate: parseFloat(row.rate) || 500,
      glr: parseFloat(row.glr) || 400,
      pwh: parseFloat(row.pwh) || 200,
      pwf: parseFloat(row.pwf) || 1800,
    };

    const res = solveNewtonRaphson(target, params, tol, 40, null);
    row.solvedField = target;

    if (target === "pwh" || target === "pwf") {
      row[target] = res.root.toFixed(1);
    } else if (target === "depth" || target === "rate" || target === "glr") {
      row[target] = Math.round(res.root).toString();
    } else if (target === "tubing") {
      row.tubing = res.root.toFixed(3);
      const nearest = getNearestApiTubing(res.root);
      row.apiNote = `API: ${nearest.od}`;
    }

    row.status = res.nonFlowing
      ? "Non-Flowing (THP=0)"
      : `Solved ${target.toUpperCase()} (${res.history.length} iters)`;
    solvedCount++;
  });

  renderBatchTable();
  const banner = document.getElementById("batchBanner");
  banner.classList.remove("hidden");
  if (skippedCount > 0) {
    banner.className =
      "text-xs px-3.5 py-2 rounded-lg border bg-amber-50 border-amber-200 text-amber-800";
    banner.textContent = `Solved ${solvedCount} row(s). ${skippedCount} row(s) skipped because more than 1 cell was blank.`;
  } else {
    banner.className =
      "text-xs px-3.5 py-2 rounded-lg border bg-emerald-50 border-emerald-200 text-emerald-800";
    banner.textContent = `Successfully solved ${solvedCount} well case(s) via Newton-Raphson.`;
  }
  showToast(`Solved ${solvedCount} batch well row(s)`);
}

function loadBatchRow(idx) {
  const row = batchRows[idx];
  const fields = ["depth", "tubing", "rate", "glr", "pwh", "pwf"];
  const blanks = fields.filter(
    (f) => row[f] === "" || isNaN(parseFloat(row[f])),
  );
  const target = row.solvedField || (blanks.length === 1 ? blanks[0] : "pwh");

  document.getElementById("targetVar").value = target;
  if (row.depth !== "") document.getElementById("depth").value = row.depth;
  if (row.tubing !== "") {
    document.getElementById("tubingID").value = row.tubing;
    syncTubingFromInput();
  }
  if (row.rate !== "") document.getElementById("rate").value = row.rate;
  if (row.glr !== "") document.getElementById("glr").value = row.glr;
  if (row.pwh !== "") document.getElementById("pwh").value = row.pwh;
  if (row.pwf !== "") document.getElementById("pwf").value = row.pwf;

  handleTargetChange();
  window.scrollTo({ top: 0, behavior: "smooth" });
  showToast(`Loaded "${row.id}" into main solver`);
}

function applyPreset(name) {
  userEditedGuess = false;
  if (name === "nind_benchmark") {
    document.getElementById("targetVar").value = "pwh";
    document.getElementById("depth").value = "8000";
    document.getElementById("tubingID").value = "1.995";
    document.getElementById("tubingPreset").value = "1.995";
    document.getElementById("rate").value = "600";
    document.getElementById("glr").value = "400";
    document.getElementById("pwf").value = "2200";
  } else if (name === "find_bhp") {
    document.getElementById("targetVar").value = "pwf";
    document.getElementById("depth").value = "7500";
    document.getElementById("tubingID").value = "1.995";
    document.getElementById("tubingPreset").value = "1.995";
    document.getElementById("rate").value = "500";
    document.getElementById("glr").value = "550";
    document.getElementById("pwh").value = "220";
  } else if (name === "find_thp") {
    document.getElementById("targetVar").value = "pwh";
    document.getElementById("depth").value = "6800";
    document.getElementById("tubingID").value = "2.441";
    document.getElementById("tubingPreset").value = "2.441";
    document.getElementById("rate").value = "850";
    document.getElementById("glr").value = "450";
    document.getElementById("pwf").value = "1950";
  } else if (name === "gas_lift") {
    document.getElementById("targetVar").value = "glr";
    document.getElementById("depth").value = "6500";
    document.getElementById("tubingID").value = "1.995";
    document.getElementById("tubingPreset").value = "1.995";
    document.getElementById("rate").value = "500";
    document.getElementById("pwh").value = "180";
    document.getElementById("pwf").value = "1650";
  } else if (name === "valve_depth") {
    document.getElementById("targetVar").value = "depth";
    document.getElementById("tubingID").value = "1.995";
    document.getElementById("tubingPreset").value = "1.995";
    document.getElementById("rate").value = "450";
    document.getElementById("glr").value = "650";
    document.getElementById("pwh").value = "160";
    document.getElementById("pwf").value = "1750";
  }
  updateNearestApiDisplay(
    parseFloat(document.getElementById("tubingID").value),
  );
  handleTargetChange();
  showToast("Loaded preset case");
}

function exportTraverseCSV() {
  if (!lastSolvedState) return;
  const s = lastSolvedState;
  let csv =
    "Station,True_Depth_ft,Flowing_Pressure_psig,Gradient_psi_per_ft,Equivalent_Depth_H_ft\n";
  const steps = 20;
  for (let i = 0; i <= steps; i++) {
    const z = (i / steps) * s.depth;
    const eqH = s.h1 + z;
    const p =
      i === 0
        ? s.pwh
        : i === steps
          ? s.pwf
          : invertEquivalentDepth(eqH, s.rate, s.glr, s.id);
    const grad = getGradient(p, s.rate, s.glr, s.id);
    csv += `${i},${Math.round(z)},${p.toFixed(2)},${grad.toFixed(4)},${Math.round(eqH)}\n`;
  }
  triggerDownload(
    csv,
    `gilbert_traverse_D${Math.round(s.depth)}_q${Math.round(s.rate)}.csv`,
  );
  showToast("Downloaded depth traverse CSV");
}

function exportBatchCSV() {
  let csv =
    "Well_Case_ID,Depth_ft,Tubing_ID_in,Rate_bbl_day,GLR_scf_bbl,THP_psig,BHP_psig,Solved_Parameter,Status\n";
  batchRows.forEach((r) => {
    csv += `"${r.id}",${r.depth},${r.tubing},${r.rate},${r.glr},${r.pwh},${r.pwf},${r.solvedField || "None"},"${r.status}"\n`;
  });
  triggerDownload(csv, "gilbert_newton_raphson_batch_cases.csv");
  showToast("Downloaded batch field data CSV");
}

function triggerDownload(content, filename) {
  const blob = new Blob([content], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.setAttribute("download", filename);
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}

let toastTimeout = null;
function showToast(msg) {
  const el = document.getElementById("toastBanner");
  document.getElementById("toastText").textContent = msg;
  el.classList.remove("translate-y-24", "opacity-0");
  el.classList.add("translate-y-0", "opacity-100");
  clearTimeout(toastTimeout);
  toastTimeout = setTimeout(() => {
    el.classList.add("translate-y-24", "opacity-0");
    el.classList.remove("translate-y-0", "opacity-100");
  }, 2400);
}

window.addEventListener("resize", drawCanvas);
window.addEventListener("DOMContentLoaded", () => {
  renderBatchTable();
  handleTargetChange();
});
