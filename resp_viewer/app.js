/**
 * CapnoView
 * - Volumetrikus kapnográfia (CO2 vs. kilégzett volumen):
 *   1. Fowler-féle anatómiai holttér (VD,Fowler) a II. fázis inflexiós pontjából
 *   2. Kilégzett átlagos CO2 (PECO2) és alveoláris átlagos CO2 (PACO2)
 *   3. Bohr-féle élettani holttér: V_DB = (PACO2 - PECO2) / PACO2 (% és ml)
 *   4. Kilégzett CO2-térfogat: VCO2 (ml/ciklus és ml/min)
 * - Utólagos kalibráció (CO2-kapnográf és pneumotachográf áramlásjel)
 * - Min-Max pixel column decimation canvas rendering (< 0.5 ms renderelési idő)
 */

class RespViewerApp {
  constructor() {
    this.rawFile = null;
    this.time_s = [];
    this.flow = [];
    this.ch4 = [];
    this.vol = [];
    this.breaths = [];

    // Alapértelmezett kalibráció (EtCO2 ~ 36 Hgmm, VT ~ 400-600 ml átlag ~ 500 ml)
    this.co2Gain = 0.231;
    this.co2Zero = 22.0;
    this.flowGain = 19.56;
    this.flowZeroOffset = 0.0;

    // Jelfeldolgozási beállítások
    this.volMode = 'linear_detrend';
    this.detectMode = 'capnogram';
    this.polarity = 1.0;
    this.baselineOffset = 0.0;
    this.viewMode = 'stacked'; // 'stacked', 'overlay', 'loop', 'volcap'

    // Interaktív nézet
    this.viewMinS = 0.0;
    this.viewMaxS = 60.0;
    this.isDragging = false;
    this.dragStartX = 0;
    this.dragStartMinS = 0;
    this.dragStartMaxS = 0;
    this.activeBreathIndex = -1;

    this.cacheDOM();
    this.bindEvents();
    this.initCanvasResize();
    this.fetchWorkspaceFiles();
  }

  cacheDOM() {
    this.fileSelect = document.getElementById('workspace-files-select');
    this.refreshBtn = document.getElementById('refresh-files-btn');
    this.localFileInput = document.getElementById('local-file-input');
    this.dropZone = document.getElementById('drop-zone');
    this.loadingOverlay = document.getElementById('loading-overlay');

    // Kalibrációs mezők
    this.calCo2GainInput = document.getElementById('cal-co2-gain');
    this.calCo2ZeroInput = document.getElementById('cal-co2-zero');
    this.calFlowGainInput = document.getElementById('cal-flow-gain');
    this.calFlowZeroInput = document.getElementById('cal-flow-zero');
    this.calResetBtn = document.getElementById('cal-reset-btn');

    // Vezérlők
    this.volMethodSelect = document.getElementById('vol-method-select');
    this.detectMethodSelect = document.getElementById('detect-method-select');
    this.polarityToggle = document.getElementById('polarity-toggle');
    this.baselineSlider = document.getElementById('baseline-slider');
    this.baselineValLabel = document.getElementById('baseline-val-label');
    this.resetZoomBtn = document.getElementById('reset-zoom-btn');
    this.modeTabs = document.querySelectorAll('.mode-tab');
    this.tooltipBar = document.getElementById('hover-tooltip-bar');
    this.visibleRangeLabel = document.getElementById('visible-range-label');
    this.volcapLegend = document.getElementById('volcap-legend');

    // Canvas elemek
    this.chartViewport = document.getElementById('chart-viewport');
    this.stackedCanvas = document.getElementById('stacked-canvas');
    this.overlayCanvas = document.getElementById('overlay-canvas');
    this.loopCanvas = document.getElementById('loop-canvas');
    this.volcapCanvas = document.getElementById('volcap-canvas');

    this.stackedCtx = this.stackedCanvas.getContext('2d');
    this.overlayCtx = this.overlayCanvas.getContext('2d');
    this.loopCtx = this.loopCanvas.getContext('2d');
    this.volcapCtx = this.volcapCanvas.getContext('2d');

    // Statisztika és táblázat
    this.statVT = document.getElementById('stat-vt');
    this.statFowler = document.getElementById('stat-fowler');
    this.statPecoPaco = document.getElementById('stat-peco-paco');
    this.statBohr = document.getElementById('stat-bohr');
    this.statVCO2 = document.getElementById('stat-vco2');
    this.statRR = document.getElementById('stat-rr');
    this.tableBody = document.getElementById('breaths-table-body');
    this.searchBreathInput = document.getElementById('breath-search-input');

    // Export
    this.exportSignalsBtn = document.getElementById('export-signals-btn');
    this.exportBreathsBtn = document.getElementById('export-breaths-btn');
  }

  bindEvents() {
    this.refreshBtn.addEventListener('click', () => this.fetchWorkspaceFiles());
    this.fileSelect.addEventListener('change', (e) => {
      if (e.target.value) this.loadFileFromWorkspace(e.target.value);
    });

    this.localFileInput.addEventListener('change', (e) => {
      if (e.target.files && e.target.files[0]) {
        this.readLocalFile(e.target.files[0]);
      }
    });

    this.dropZone.addEventListener('dragover', (e) => {
      e.preventDefault();
      this.dropZone.classList.add('drag-over');
    });
    this.dropZone.addEventListener('dragleave', () => {
      this.dropZone.classList.remove('drag-over');
    });
    this.dropZone.addEventListener('drop', (e) => {
      e.preventDefault();
      this.dropZone.classList.remove('drag-over');
      if (e.dataTransfer.files && e.dataTransfer.files[0]) {
        this.readLocalFile(e.dataTransfer.files[0]);
      }
    });

    // Kalibrációs események
    this.calCo2GainInput.addEventListener('input', (e) => {
      this.co2Gain = parseFloat(e.target.value) || 0.231;
      this.reprocessSignals();
    });
    this.calCo2ZeroInput.addEventListener('input', (e) => {
      this.co2Zero = parseFloat(e.target.value) || 22.0;
      this.reprocessSignals();
    });
    this.calFlowGainInput.addEventListener('input', (e) => {
      this.flowGain = parseFloat(e.target.value) || 19.56;
      this.reprocessSignals();
    });
    this.calFlowZeroInput.addEventListener('input', (e) => {
      this.flowZeroOffset = parseFloat(e.target.value) || 0.0;
      this.reprocessSignals();
    });

    this.calResetBtn.addEventListener('click', () => {
      this.co2Gain = 0.231;
      this.co2Zero = 22.0;
      this.flowGain = 19.56;
      this.flowZeroOffset = 0.0;
      this.calCo2GainInput.value = "0.231";
      this.calCo2ZeroInput.value = "22.0";
      this.calFlowGainInput.value = "19.56";
      this.calFlowZeroInput.value = "0.0";
      this.reprocessSignals();
    });

    this.volMethodSelect.addEventListener('change', (e) => {
      this.volMode = e.target.value;
      this.reprocessSignals();
    });

    this.detectMethodSelect.addEventListener('change', (e) => {
      this.detectMode = e.target.value;
      this.reprocessSignals();
    });

    this.polarityToggle.addEventListener('change', (e) => {
      this.polarity = e.target.checked ? -1.0 : 1.0;
      this.reprocessSignals();
    });

    this.baselineSlider.addEventListener('input', (e) => {
      this.baselineOffset = parseFloat(e.target.value);
      this.baselineValLabel.textContent = this.baselineOffset.toFixed(1);
      this.reprocessSignals();
    });

    this.resetZoomBtn.addEventListener('click', () => this.resetZoom());

    this.modeTabs.forEach(tab => {
      tab.addEventListener('click', (e) => {
        this.modeTabs.forEach(t => t.classList.remove('active'));
        e.target.classList.add('active');
        this.setViewMode(e.target.dataset.mode);
      });
    });

    [this.stackedCanvas, this.overlayCanvas, this.loopCanvas, this.volcapCanvas].forEach(canvas => {
      canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e, canvas));
      canvas.addEventListener('mouseleave', () => {
        this.tooltipBar.textContent = "Mozgass egeret a grafikonok felett a pontos idő-, áramlás (ml/s)-, CO2 (Hgmm)- és volumenértékek (ml) kiolvasásához...";
        this.renderCharts();
      });
      canvas.addEventListener('mousedown', (e) => this.handleMouseDown(e));
      canvas.addEventListener('mouseup', () => this.handleMouseUp());
      canvas.addEventListener('wheel', (e) => this.handleWheel(e), { passive: false });
      canvas.addEventListener('dblclick', () => this.resetZoom());
    });

    this.searchBreathInput.addEventListener('input', () => this.renderBreathsTable());
    this.exportSignalsBtn.addEventListener('click', () => this.exportSignalsCSV());
    this.exportBreathsBtn.addEventListener('click', () => this.exportBreathsCSV());
  }

  showLoading(text) {
    document.getElementById('loading-text').textContent = text;
    this.loadingOverlay.classList.remove('hidden');
  }

  hideLoading() {
    this.loadingOverlay.classList.add('hidden');
  }

  async fetchWorkspaceFiles() {
    try {
      const resp = await fetch('/api/files');
      if (!resp.ok) return;
      const data = await resp.json();
      this.fileSelect.innerHTML = '<option value="">-- Válassz fájlt a mappából --</option>';
      data.files.forEach(f => {
        const opt = document.createElement('option');
        opt.value = f.name;
        opt.textContent = `${f.name} (${f.type.toUpperCase()}, ${f.size_kb} KB)`;
        if (f.name === 'pc002_002.inp') {
          opt.selected = true;
        }
        this.fileSelect.appendChild(opt);
      });

      if (this.fileSelect.value) {
        this.loadFileFromWorkspace(this.fileSelect.value);
      }
    } catch (e) {
      console.log('Statikus mód vagy helyi fájlelérés');
    }
  }

  async loadFileFromWorkspace(filename) {
    this.showLoading(`Betöltés a szerverről: ${filename}...`);
    try {
      const resp = await fetch(`/api/load?file=${encodeURIComponent(filename)}`);
      if (!resp.ok) throw new Error('Hiba a fájl betöltésekor');
      const data = await resp.json();
      this.processLoadedData(data);
    } catch (err) {
      alert(`Hiba történt: ${err.message}`);
    } finally {
      this.hideLoading();
    }
  }

  async readLocalFile(file) {
    this.showLoading(`Helyi fájl olvasása: ${file.name}...`);
    try {
      const ext = file.name.split('.').pop().toLowerCase();
      if (ext === 'inp') {
        const buffer = await file.arrayBuffer();
        const data = this.parseInpBuffer(buffer, file.name);
        this.processLoadedData(data);
      } else if (ext === 'csv') {
        const text = await file.text();
        const data = this.parseCsvText(text, file.name);
        this.processLoadedData(data);
      }
    } catch (err) {
      alert(`Hiba a fájl feldolgozása közben: ${err.message}`);
    } finally {
      this.hideLoading();
    }
  }

  parseInpBuffer(buffer, filename) {
    const view = new DataView(buffer);
    const sampleRateHz = view.getFloat32(0, true);
    const sampleCount = view.getUint16(4, true);
    const channelCount = view.getUint16(6, true);
    const durationS = view.getFloat32(8, true);

    let pos = 42;
    const readPascal = () => {
      const len = view.getUint16(pos, true);
      pos += 2;
      let str = '';
      for (let i = 0; i < len; i++) {
        str += String.fromCharCode(view.getUint8(pos + i));
      }
      pos += len;
      return str;
    };

    const recordedAt = readPascal();
    const subjectId = readPascal();
    pos += 2;
    const sgnPath = readPascal();

    const int16Array = new Int16Array(buffer, pos);
    const time_s = [];
    const channels = {};
    for (let c = 0; c < channelCount; c++) {
      channels[`raw_ch${c + 1}`] = [];
    }

    let idx = 0;
    for (let s = 0; s < sampleCount; s++) {
      time_s.push(s / sampleRateHz);
      for (let c = 0; c < channelCount; c++) {
        channels[`raw_ch${c + 1}`].push(int16Array[idx++]);
      }
    }

    return {
      filename,
      file_type: 'inp',
      header: {
        sample_rate_hz: sampleRateHz,
        sample_count: sampleCount,
        channel_count: channelCount,
        duration_s: durationS,
        recorded_at: recordedAt,
        subject_id: subjectId
      },
      time_s,
      channels
    };
  }

  parseCsvText(text, filename) {
    const lines = text.trim().split(/\r?\n/);
    const headers = lines[0].split(',').map(h => h.trim());
    const dataRows = lines.slice(1).map(l => l.split(',').map(Number));

    const time_s = dataRows.map(r => r[0]);
    const channels = {};

    headers.slice(1).forEach((h, colIdx) => {
      channels[h] = dataRows.map(r => r[colIdx + 1]);
      if (h.toLowerCase().includes('ch1') || h.toLowerCase().includes('respiratory')) {
        channels['raw_ch1'] = channels[h];
      }
      if (h.toLowerCase().includes('ch4') || h.toLowerCase().includes('capnogram') || h.toLowerCase().includes('co2')) {
        channels['raw_ch4'] = channels[h];
      }
    });

    const sampleRateHz = time_s.length > 1 ? 1.0 / (time_s[1] - time_s[0]) : 256.0;

    return {
      filename,
      file_type: 'csv',
      header: {
        sample_rate_hz: sampleRateHz,
        sample_count: time_s.length,
        channel_count: headers.length - 1,
        duration_s: time_s[time_s.length - 1] || 0,
        recorded_at: 'CSV adatsor',
        subject_id: filename.replace(/\.[^/.]+$/, "")
      },
      time_s,
      channels
    };
  }

  processLoadedData(data) {
    this.rawFile = data;
    this.time_s = data.time_s;
    this.ch1Raw = data.channels['raw_ch1'] || [];
    this.ch4Raw = data.channels['raw_ch4'] || null;

    this.viewMinS = this.time_s[0] || 0.0;
    this.viewMaxS = this.time_s[this.time_s.length - 1] || 60.0;

    // Automatic polarity inversion for pc1 directory files
    const filename = data.filename || '';
    if (data.default_polarity !== undefined) {
      this.polarity = data.default_polarity;
    } else if (filename.includes('pc1/') || filename.startsWith('pc1')) {
      this.polarity = -1.0;
    } else {
      this.polarity = 1.0;
    }
    if (this.polarityToggle) {
      this.polarityToggle.checked = (this.polarity === -1.0);
    }

    this.updateMetadataUI(data);

    this.exportSignalsBtn.disabled = false;
    this.exportBreathsBtn.disabled = false;

    this.reprocessSignals();
  }

  updateMetadataUI(data) {
    document.getElementById('file-type-badge').textContent = data.file_type.toUpperCase();
    document.getElementById('file-type-badge').className = `badge ${data.file_type}`;
    document.getElementById('meta-subject').textContent = data.header.subject_id || '-';
    document.getElementById('meta-date').textContent = data.header.recorded_at || '-';
    document.getElementById('meta-fs').textContent = `${Math.round(data.header.sample_rate_hz)} Hz`;
    document.getElementById('meta-duration').textContent = `${data.header.duration_s.toFixed(1)} s`;
    document.getElementById('meta-channels').textContent = `${data.header.channel_count} csatorna`;
  }

  reprocessSignals() {
    if (!this.rawFile || !this.time_s.length) return;

    const n = this.time_s.length;
    const dt = this.time_s[1] - this.time_s[0];

    // 1. Kalibrált áramlásjel ml/s
    // Auto-balanced mean baseline: mathematically forces Vinsp == Vexp (zero net volume drift)
    let sumCh1 = 0;
    for (let i = 0; i < n; i++) sumCh1 += this.ch1Raw[i];
    const meanCh1 = sumCh1 / n;
    const baseline = meanCh1 + this.baselineOffset + this.flowZeroOffset;

    this.flow = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      this.flow[i] = -(this.ch1Raw[i] - baseline) * this.polarity * this.flowGain;
    }

    // 2. Kalibrált kapnogram Hgmm
    this.ch4 = new Float32Array(n);
    if (this.ch4Raw) {
      for (let i = 0; i < n; i++) {
        this.ch4[i] = (this.ch4Raw[i] - this.co2Zero) * this.co2Gain;
      }
    }

    // 3. Adaptív ciklusérzékelés
    this.breathStarts = this.detectBreathBoundaries(this.ch4, this.flow, this.time_s, this.detectMode);

    // 4. Számított kalibrált volumenjel ml
    this.vol = new Float32Array(n);
    const volRaw = new Float32Array(n);
    let accum = 0;
    for (let i = 0; i < n; i++) {
      accum += this.flow[i] * dt;
      volRaw[i] = accum;
    }

    if (this.volMode === 'linear_detrend') {
      for (let k = 0; k < this.breathStarts.length - 1; k++) {
        const idx0 = this.breathStarts[k];
        const idx1 = this.breathStarts[k + 1];
        const count = idx1 - idx0 + 1;
        const drift = volRaw[idx1] - volRaw[idx0];
        const startVal = volRaw[idx0];

        for (let j = 0; j < count; j++) {
          const frac = j / (count - 1);
          const correction = startVal + drift * frac;
          const val = volRaw[idx0 + j] - correction;
          if (k === 0 || j > 0) {
            this.vol[idx0 + j] = val;
          }
        }
      }
    } else if (this.volMode === 'envelope') {
      for (let i = 0; i < n; i++) {
        let idxPrev = this.breathStarts[0];
        let idxNext = this.breathStarts[this.breathStarts.length - 1];
        for (let k = 0; k < this.breathStarts.length - 1; k++) {
          if (i >= this.breathStarts[k] && i <= this.breathStarts[k + 1]) {
            idxPrev = this.breathStarts[k];
            idxNext = this.breathStarts[k + 1];
            break;
          }
        }
        const frac = idxNext > idxPrev ? (i - idxPrev) / (idxNext - idxPrev) : 0;
        const envelope = volRaw[idxPrev] + (volRaw[idxNext] - volRaw[idxPrev]) * frac;
        this.vol[i] = volRaw[i] - envelope;
      }
    } else if (this.volMode === 'highpass') {
      let alpha = 0.995;
      let prevVal = 0;
      for (let i = 0; i < n; i++) {
        const currVal = alpha * prevVal + this.flow[i] * dt;
        this.vol[i] = currVal;
        prevVal = currVal;
      }
    } else {
      for (let i = 0; i < n; i++) {
        this.vol[i] = volRaw[i];
      }
    }

    // 5. Ciklusok és volumetrikus kapnográfia paramétereinek számítása
    this.breaths = this.calculateBreathStatistics(this.breathStarts, this.time_s, this.vol, this.ch4);

    this.updateStatsCards();
    this.renderBreathsTable();
    this.renderCharts();
  }

  detectBreathBoundaries(ch4, flow, time_s, mode) {
    const n = time_s.length;
    if (!n || !ch4 || ch4.length !== n) return [0, n - 1];

    const dt = time_s[1] - time_s[0];
    let minCo2 = Infinity, maxCo2 = -Infinity;
    for (let i = 0; i < n; i++) {
      if (ch4[i] < minCo2) minCo2 = ch4[i];
      if (ch4[i] > maxCo2) maxCo2 = ch4[i];
    }
    const co2Range = maxCo2 - minCo2;
    if (co2Range < 5.0) return [0, n - 1];

    // Step 1: Coarse CO2 expiration boundaries
    const threshHigh = minCo2 + co2Range * 0.25;
    const threshLow = minCo2 + co2Range * 0.15;
    const coarseStarts = [];
    const coarseEnds = [];
    let inExp = false;

    for (let i = 1; i < n; i++) {
      if (!inExp && ch4[i] > threshHigh) {
        if (coarseStarts.length === 0 || (time_s[i] - time_s[coarseStarts[coarseStarts.length - 1]]) > 0.8) {
          coarseStarts.push(i);
          inExp = true;
        }
      } else if (inExp && ch4[i] < threshLow) {
        coarseEnds.push(i);
        inExp = false;
      }
    }

    const nCoarse = Math.min(coarseStarts.length, coarseEnds.length);
    if (nCoarse === 0) return [0, n - 1];

    // Step 2: Fine-tuning via calibrated flow zero-crossings (+/- 0.5s window)
    const winSamples = Math.floor(0.5 / dt);
    const starts = [];

    for (let k = 0; k < nCoarse; k++) {
      const cs = coarseStarts[k];
      const ce = coarseEnds[k];

      const w0 = Math.max(1, cs - winSamples);
      const w1 = Math.min(n - 1, cs + winSamples);
      let fs = cs;
      for (let i = w0; i < w1; i++) {
        if (flow[i - 1] <= 0 && flow[i] > 0) {
          fs = i;
          break;
        }
      }

      const w2 = Math.max(1, ce - winSamples);
      const w3 = Math.min(n - 1, ce + winSamples);
      let fe = ce;
      for (let i = w2; i < w3; i++) {
        if (flow[i - 1] >= 0 && flow[i] < 0) {
          fe = i;
          break;
        }
      }

      if (fe > fs + 10) {
        starts.push(fs);
      }
    }

    if (starts.length === 0) return [0, n - 1];
    if (starts[starts.length - 1] !== n - 1) starts.push(n - 1);

    return starts;
  }

  calculateBreathStatistics(starts, time_s, vol, ch4) {
    const list = [];
    for (let k = 0; k < starts.length - 1; k++) {
      const i0 = starts[k];
      const i1 = starts[k + 1];

      let maxVol = -Infinity;
      let peakIdx = i0;
      for (let j = i0; j <= i1; j++) {
        if (vol[j] > maxVol) {
          maxVol = vol[j];
          peakIdx = j;
        }
      }

      const ti = time_s[peakIdx] - time_s[i0];
      const te = time_s[i1] - time_s[peakIdx];
      const totTime = ti + te;
      const rr = totTime > 0 ? 60.0 / totTime : 0.0;
      const ieRatioVal = ti > 0 ? (te / ti) : 0;

      let etCo2 = 0;
      let fowler_vd = 0;
      let peco2 = 0;
      let paco2 = 0;
      let bohr_ratio = 0;
      let bohr_vd = 0;
      let vco2_ml = 0;
      const volcap_curve = [];

      const vt = maxVol;

      if (ch4 && ch4.length > i1 && vt > 50 && (i1 - peakIdx) > 5) {
        const co2_exp = [];
        const v_exh = [];
        for (let j = peakIdx; j <= i1; j++) {
          co2_exp.push(ch4[j]);
          v_exh.push(vt - vol[j]);
          if (ch4[j] > etCo2) etCo2 = ch4[j];
        }

        const maxVInt = Math.floor(vt);
        if (maxVInt > 10) {
          const indices = v_exh.map((v, idx) => ({ v, idx })).sort((a, b) => a.v - b.v);
          const v_sorted = [];
          const co2_sorted = [];
          for (let idx = 0; idx < indices.length; idx++) {
            const val = indices[idx].v;
            if (v_sorted.length === 0 || Math.abs(val - v_sorted[v_sorted.length - 1]) > 0.01) {
              v_sorted.push(val);
              co2_sorted.push(co2_exp[indices[idx].idx]);
            }
          }

          const v_grid = [];
          const co2_grid = [];
          for (let g = 0; g <= maxVInt; g++) {
            v_grid.push(g);
            co2_grid.push(this.linearInterp(g, v_sorted, co2_sorted));
          }

          const stepDs = Math.max(1, Math.floor(v_grid.length / 80));
          for (let g = 0; g < v_grid.length; g += stepDs) {
            volcap_curve.push({ v: v_grid[g], co2: co2_grid[g] });
          }

          let sumCo2 = 0;
          for (let g = 0; g < co2_grid.length; g++) sumCo2 += co2_grid[g];
          peco2 = sumCo2 / co2_grid.length;

          const dco2 = [];
          for (let g = 0; g < co2_grid.length; g++) {
            const g0 = Math.max(0, g - 5);
            const g1 = Math.min(co2_grid.length - 1, g + 5);
            dco2.push((co2_grid[g1] - co2_grid[g0]) / (g1 - g0 || 1));
          }
          const iMin = Math.floor(0.15 * v_grid.length);
          const iMax = Math.floor(0.65 * v_grid.length);
          let maxSlope = -1, idxInflect = iMin;
          for (let g = iMin; g < iMax && g < dco2.length; g++) {
            if (dco2[g] > maxSlope) {
              maxSlope = dco2[g];
              idxInflect = g;
            }
          }
          fowler_vd = v_grid[idxInflect] || (0.3 * vt);

          const iP3 = Math.min(Math.floor(fowler_vd + 0.15 * vt), Math.floor(0.65 * v_grid.length));
          let sumP3 = 0, countP3 = 0;
          for (let g = iP3; g < co2_grid.length; g++) {
            sumP3 += co2_grid[g];
            countP3++;
          }
          paco2 = countP3 > 0 ? (sumP3 / countP3) : etCo2;

          bohr_ratio = paco2 > 0 ? ((paco2 - peco2) / paco2) : 0;
          bohr_vd = bohr_ratio * vt;

          vco2_ml = (peco2 / 760.0) * vt;
        }
      }

      list.push({
        index: k + 1,
        startIdx: i0,
        peakIdx: peakIdx,
        endIdx: i1,
        start_s: time_s[i0],
        peak_s: time_s[peakIdx],
        end_s: time_s[i1],
        vt: maxVol,
        ti: ti,
        te: te,
        ieRatioVal: ieRatioVal,
        rr: rr,
        et_co2: etCo2,
        fowler_vd: fowler_vd,
        peco2: peco2,
        paco2: paco2,
        bohr_ratio: bohr_ratio,
        bohr_vd: bohr_vd,
        vco2_ml: vco2_ml,
        volcap_curve: volcap_curve
      });
    }
    return list;
  }

  linearInterp(x, xs, ys) {
    if (!xs.length) return 0;
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    for (let i = 0; i < xs.length - 1; i++) {
      if (x >= xs[i] && x <= xs[i + 1]) {
        const frac = (x - xs[i]) / (xs[i + 1] - xs[i] || 1);
        return ys[i] + (ys[i + 1] - ys[i]) * frac;
      }
    }
    return ys[0];
  }

  updateStatsCards() {
    document.getElementById('meta-breaths-count').textContent = `${this.breaths.length} ciklus`;

    if (!this.breaths.length) {
      this.statVT.innerHTML = `0.0 <small>ml</small>`;
      this.statFowler.innerHTML = `0.0 <small>ml</small>`;
      this.statPecoPaco.innerHTML = `0.0 / 0.0 <small>Hgmm</small>`;
      this.statBohr.innerHTML = `0.0 <small>%</small>`;
      this.statVCO2.innerHTML = `0.0 <small>ml</small>`;
      this.statRR.innerHTML = `0.0 <small>bpm</small>`;
      return;
    }

    const meanVT = this.breaths.reduce((s, b) => s + b.vt, 0) / this.breaths.length;
    const meanFowler = this.breaths.reduce((s, b) => s + b.fowler_vd, 0) / this.breaths.length;
    const meanPECO2 = this.breaths.reduce((s, b) => s + b.peco2, 0) / this.breaths.length;
    const meanPACO2 = this.breaths.reduce((s, b) => s + b.paco2, 0) / this.breaths.length;
    const meanBohrR = this.breaths.reduce((s, b) => s + b.bohr_ratio, 0) / this.breaths.length;
    const meanVCO2 = this.breaths.reduce((s, b) => s + b.vco2_ml, 0) / this.breaths.length;
    const meanRR = this.breaths.reduce((s, b) => s + b.rr, 0) / this.breaths.length;

    this.statVT.innerHTML = `${meanVT.toFixed(1)} <small>ml</small>`;
    this.statFowler.innerHTML = `${meanFowler.toFixed(1)} <small>ml</small>`;
    this.statPecoPaco.innerHTML = `${meanPECO2.toFixed(1)} / ${meanPACO2.toFixed(1)} <small>Hgmm</small>`;
    this.statBohr.innerHTML = `${(meanBohrR * 100).toFixed(1)} <small>%</small>`;
    this.statVCO2.innerHTML = `${meanVCO2.toFixed(1)} <small>ml</small>`;
    this.statRR.innerHTML = `${meanRR.toFixed(1)} <small>bpm</small>`;
  }

  renderBreathsTable() {
    this.tableBody.innerHTML = '';
    const query = this.searchBreathInput.value.trim().toLowerCase();

    const filtered = this.breaths.filter(b => {
      if (!query) return true;
      return b.index.toString().includes(query);
    });

    if (!filtered.length) {
      this.tableBody.innerHTML = '<tr class="empty-row"><td colspan="14">Nincs a szűrésnek megfelelő légzési ciklus.</td></tr>';
      return;
    }

    const frag = document.createDocumentFragment();
    filtered.forEach(b => {
      const tr = document.createElement('tr');
      if (b.index === this.activeBreathIndex) {
        tr.classList.add('active-row');
      }
      tr.innerHTML = `
        <td><strong>#${b.index}</strong></td>
        <td><strong style="color:var(--color-vol)">${b.vt.toFixed(1)}</strong></td>
        <td><span style="color:var(--color-co2)">${b.et_co2.toFixed(1)}</span></td>
        <td><strong style="color:#00f0ff">${b.fowler_vd.toFixed(1)}</strong></td>
        <td>${b.peco2.toFixed(1)}</td>
        <td><strong>${b.paco2.toFixed(1)}</strong></td>
        <td><span style="color:#f43f5e">${(b.bohr_ratio * 100).toFixed(1)}%</span></td>
        <td>${b.bohr_vd.toFixed(1)}</td>
        <td><strong style="color:#10b981">${b.vco2_ml.toFixed(1)}</strong></td>
        <td>${b.start_s.toFixed(2)}</td>
        <td>${b.end_s.toFixed(2)}</td>
        <td>${b.ti.toFixed(2)}</td>
        <td>${b.te.toFixed(2)}</td>
        <td>${b.rr.toFixed(1)}</td>
      `;

      tr.addEventListener('click', () => {
        this.activeBreathIndex = b.index;
        this.zoomToBreath(b);
        if (this.viewMode !== 'volcap') {
          this.setViewMode('volcap');
          this.modeTabs.forEach(t => t.classList.toggle('active', t.dataset.mode === 'volcap'));
        }
        this.renderBreathsTable();
      });

      frag.appendChild(tr);
    });
    this.tableBody.appendChild(frag);
  }

  zoomToBreath(b) {
    const margin = (b.end_s - b.start_s) * 0.3;
    this.viewMinS = Math.max(0, b.start_s - margin);
    this.viewMaxS = Math.min(this.time_s[this.time_s.length - 1], b.end_s + margin);
    this.renderCharts();
  }

  resetZoom() {
    if (!this.time_s.length) return;
    this.viewMinS = this.time_s[0];
    this.viewMaxS = this.time_s[this.time_s.length - 1];
    this.activeBreathIndex = -1;
    this.renderBreathsTable();
    this.renderCharts();
  }

  setViewMode(mode) {
    this.viewMode = mode;
    this.stackedCanvas.classList.toggle('hidden', mode !== 'stacked');
    this.overlayCanvas.classList.toggle('hidden', mode !== 'overlay');
    this.loopCanvas.classList.toggle('hidden', mode !== 'loop');
    this.volcapCanvas.classList.toggle('hidden', mode !== 'volcap');
    this.volcapLegend.classList.toggle('hidden', mode !== 'volcap');
    this.renderCharts();
  }

  initCanvasResize() {
    const resizeObserver = new ResizeObserver(() => {
      this.resizeCanvases();
      this.renderCharts();
    });
    resizeObserver.observe(this.chartViewport);
  }

  resizeCanvases() {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.chartViewport.getBoundingClientRect();
    const width = rect.width;
    const height = rect.height;

    [this.stackedCanvas, this.overlayCanvas, this.loopCanvas, this.volcapCanvas].forEach(canvas => {
      canvas.width = width * dpr;
      canvas.height = height * dpr;
      canvas.style.width = `${width}px`;
      canvas.style.height = `${height}px`;
      const ctx = canvas.getContext('2d');
      ctx.resetTransform();
      ctx.scale(dpr, dpr);
    });
  }

  renderCharts() {
    if (!this.time_s.length) return;
    this.visibleRangeLabel.textContent = `${this.viewMinS.toFixed(2)}s — ${this.viewMaxS.toFixed(2)}s`;

    if (this.viewMode === 'stacked') {
      this.renderStackedChart();
    } else if (this.viewMode === 'overlay') {
      this.renderOverlayChart();
    } else if (this.viewMode === 'loop') {
      this.renderLoopChart();
    } else if (this.viewMode === 'volcap') {
      this.renderVolcapChart();
    }
  }

  renderStackedChart() {
    const ctx = this.stackedCtx;
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const h = rect.height;

    ctx.clearRect(0, 0, w, h);
    const trackH = h / 3;

    this.drawPhaseBands(ctx, w, h);

    // Track 1: Airflow (ml/s)
    this.drawSignalTrack(ctx, this.flow, 0, trackH, 'var(--color-flow)', 'Spontán légzési áramlás [1. csatorna] — ml/s', true, true);
    // Track 2: Capnogram (Hgmm)
    this.drawSignalTrack(ctx, this.ch4, trackH, trackH, 'var(--color-co2)', 'Mainstream időkapnogram [4. csatorna] — Hgmm', false, false);
    // Track 3: Volume (ml)
    this.drawSignalTrack(ctx, this.vol, trackH * 2, trackH, 'var(--color-vol)', 'Számított légzési volumenjel (minden kilégzés végén 0) — ml', true, false, true);

    ctx.strokeStyle = 'rgba(255, 255, 255, 0.12)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, trackH); ctx.lineTo(w, trackH);
    ctx.moveTo(0, trackH * 2); ctx.lineTo(w, trackH * 2);
    ctx.stroke();
  }

  renderOverlayChart() {
    const ctx = this.overlayCtx;
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const h = rect.height;

    ctx.clearRect(0, 0, w, h);
    this.drawPhaseBands(ctx, w, h);

    this.drawOverlaySignal(ctx, this.flow, 0, h, 'var(--color-flow)', 0.5);
    this.drawOverlaySignal(ctx, this.ch4, 0, h, 'var(--color-co2)', 0.4);
    this.drawOverlaySignal(ctx, this.vol, 0, h, 'var(--color-vol)', 1.0, true);
  }

  renderLoopChart() {
    const ctx = this.loopCtx;
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const h = rect.height;

    ctx.clearRect(0, 0, w, h);

    ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, h / 2); ctx.lineTo(w, h / 2);
    ctx.moveTo(w / 10, 0); ctx.lineTo(w / 10, h);
    ctx.stroke();

    ctx.fillStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.font = '11px Inter';
    ctx.fillText('Légvételi volumen VT (ml) ->', w / 10 + 10, h / 2 - 10);
    ctx.fillText('Belégzési áramlás + (ml/s)', w / 10 + 10, 25);
    ctx.fillText('Kilégzési áramlás - (ml/s)', w / 10 + 10, h - 15);

    const range = this.getVisibleIndices();
    let maxV = 1, maxF = 1;
    for (let i = range.start; i <= range.end; i++) {
      if (this.vol[i] > maxV) maxV = this.vol[i];
      if (Math.abs(this.flow[i]) > maxF) maxF = Math.abs(this.flow[i]);
    }

    ctx.strokeStyle = 'var(--color-flow)';
    ctx.lineWidth = 1.8;
    ctx.beginPath();

    for (let i = range.start; i <= range.end; i++) {
      const x = w / 10 + (this.vol[i] / maxV) * (w * 0.85);
      const y = h / 2 - (this.flow[i] / maxF) * (h * 0.42);
      if (i === range.start) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  renderVolcapChart() {
    const ctx = this.volcapCtx;
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const h = rect.height;

    ctx.clearRect(0, 0, w, h);

    const padLeft = 65, padRight = 30, padTop = 40, padBottom = 45;
    const drawW = w - padLeft - padRight;
    const drawH = h - padTop - padBottom;

    ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padLeft, padTop);
    ctx.lineTo(padLeft, h - padBottom);
    ctx.lineTo(w - padRight, h - padBottom);
    ctx.stroke();

    let maxV = 600, maxCO2 = 42;
    const visibleBreaths = this.breaths.filter(b => b.end_s >= this.viewMinS && b.start_s <= this.viewMaxS);
    if (visibleBreaths.length) {
      maxV = Math.max(500, Math.max(...visibleBreaths.map(b => b.vt)) * 1.08);
      maxCO2 = Math.max(40, Math.max(...visibleBreaths.map(b => b.paco2 || b.et_co2)) * 1.15);
    }

    ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
    ctx.font = '11px Inter';
    for (let v = 0; v <= maxV; v += 100) {
      const x = padLeft + (v / maxV) * drawW;
      ctx.fillText(`${v}`, x - 10, h - padBottom + 18);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
      ctx.beginPath();
      ctx.moveTo(x, padTop);
      ctx.lineTo(x, h - padBottom);
      ctx.stroke();
    }
    for (let c = 0; c <= maxCO2; c += 10) {
      const y = h - padBottom - (c / maxCO2) * drawH;
      ctx.fillText(`${c}`, padLeft - 28, y + 4);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
      ctx.beginPath();
      ctx.moveTo(padLeft, y);
      ctx.lineTo(w - padRight, y);
      ctx.stroke();
    }

    ctx.fillStyle = 'rgba(255, 255, 255, 0.85)';
    ctx.font = '600 12px Inter';
    ctx.fillText('Kilégzett térfogat (V_exh) — ml', w / 2 - 80, h - 12);
    ctx.save();
    ctx.translate(16, h / 2 + 60);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('CO2 parciális nyomás (P_CO2) — Hgmm', 0, 0);
    ctx.restore();

    visibleBreaths.forEach(b => {
      if (!b.volcap_curve || !b.volcap_curve.length) return;
      const isActive = (b.index === this.activeBreathIndex) || (this.activeBreathIndex === -1 && b === visibleBreaths[0]);

      ctx.strokeStyle = isActive ? '#ff9f1c' : 'rgba(255, 159, 28, 0.25)';
      ctx.lineWidth = isActive ? 2.8 : 1.2;
      ctx.beginPath();
      b.volcap_curve.forEach((pt, idx) => {
        const x = padLeft + (pt.v / maxV) * drawW;
        const y = h - padBottom - (pt.co2 / maxCO2) * drawH;
        if (idx === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      });
      ctx.stroke();

      if (isActive) {
        const xFowler = padLeft + (b.fowler_vd / maxV) * drawW;
        ctx.strokeStyle = '#00f0ff';
        ctx.lineWidth = 1.8;
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(xFowler, padTop);
        ctx.lineTo(xFowler, h - padBottom);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#00f0ff';
        ctx.font = '600 11px Inter';
        ctx.fillText(`1. VD Fowler: ${b.fowler_vd.toFixed(1)} ml`, xFowler + 6, padTop + 20);

        const yPECO2 = h - padBottom - (b.peco2 / maxCO2) * drawH;
        const yPACO2 = h - padBottom - (b.paco2 / maxCO2) * drawH;

        ctx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(padLeft, yPECO2);
        ctx.lineTo(w - padRight, yPECO2);
        ctx.moveTo(padLeft, yPACO2);
        ctx.lineTo(w - padRight, yPACO2);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#f8fafc';
        ctx.fillText(`2. PECO2: ${b.peco2.toFixed(1)} Hgmm`, padLeft + 10, yPECO2 - 5);
        ctx.fillText(`2. PACO2: ${b.paco2.toFixed(1)} Hgmm (alveoláris)`, padLeft + 10, yPACO2 - 5);

        const badgeW = 260, badgeH = 58;
        const bx = w - padRight - badgeW - 15, by = padTop + 15;
        ctx.fillStyle = 'rgba(16, 24, 38, 0.88)';
        ctx.strokeStyle = 'rgba(0, 240, 255, 0.3)';
        ctx.lineWidth = 1;
        ctx.fillRect(bx, by, badgeW, badgeH);
        ctx.strokeRect(bx, by, badgeW, badgeH);

        ctx.fillStyle = '#00f0ff';
        ctx.font = '700 12px Inter';
        ctx.fillText(`Ciklus #${b.index} — Volumetrikus kapnogram`, bx + 10, by + 18);
        ctx.fillStyle = '#f8fafc';
        ctx.font = '11px Inter';
        ctx.fillText(`3. Bohr élettani holttér: ${(b.bohr_ratio * 100).toFixed(1)}% (${b.bohr_vd.toFixed(1)} ml)`, bx + 10, by + 34);
        ctx.fillStyle = '#10b981';
        ctx.fillText(`4. Kilégzett CO2 (VCO2): ${b.vco2_ml.toFixed(1)} ml/ciklus`, bx + 10, by + 50);
      }
    });
  }

  drawPhaseBands(ctx, w, h) {
    const duration = this.viewMaxS - this.viewMinS;
    for (let k = 0; k < this.breaths.length; k++) {
      const b = this.breaths[k];
      if (b.end_s < this.viewMinS || b.start_s > this.viewMaxS) continue;

      const x0 = ((b.start_s - this.viewMinS) / duration) * w;
      const xPeak = ((b.peak_s - this.viewMinS) / duration) * w;
      const xEnd = ((b.end_s - this.viewMinS) / duration) * w;

      ctx.fillStyle = 'rgba(0, 240, 255, 0.04)';
      ctx.fillRect(x0, 0, xPeak - x0, h);

      ctx.fillStyle = 'rgba(255, 159, 28, 0.04)';
      ctx.fillRect(xPeak, 0, xEnd - xPeak, h);

      if (b.index === this.activeBreathIndex) {
        ctx.strokeStyle = 'rgba(0, 240, 255, 0.5)';
        ctx.lineWidth = 2;
        ctx.strokeRect(x0, 0, xEnd - x0, h);
      }
    }
  }

  drawSignalTrack(ctx, arr, yOffset, trackH, color, label, drawZeroLine = false, isFlow = false, isVol = false) {
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const range = this.getVisibleIndices();

    let minVal = Infinity, maxVal = -Infinity;
    for (let i = range.start; i <= range.end; i++) {
      if (arr[i] < minVal) minVal = arr[i];
      if (arr[i] > maxVal) maxVal = arr[i];
    }
    if (isVol) {
      minVal = Math.min(minVal, -20);
      maxVal = Math.max(maxVal, 650);
    }
    const valRange = (maxVal - minVal) || 1;
    const paddingH = trackH * 0.15;
    const drawH = trackH - paddingH * 2;

    if (drawZeroLine) {
      const zeroY = yOffset + trackH - paddingH - ((0 - minVal) / valRange) * drawH;
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(0, zeroY);
      ctx.lineTo(w, zeroY);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    ctx.strokeStyle = color;
    ctx.lineWidth = 1.6;
    ctx.beginPath();

    const duration = this.viewMaxS - this.viewMinS;
    const numPixels = Math.ceil(w);
    for (let px = 0; px < numPixels; px++) {
      const t0 = this.viewMinS + (px / numPixels) * duration;
      const t1 = this.viewMinS + ((px + 1) / numPixels) * duration;

      const idx0 = this.findTimeIndex(t0);
      const idx1 = Math.max(idx0, this.findTimeIndex(t1));

      let sMin = arr[idx0], sMax = arr[idx0];
      for (let j = idx0; j <= idx1 && j < arr.length; j++) {
        if (arr[j] < sMin) sMin = arr[j];
        if (arr[j] > sMax) sMax = arr[j];
      }

      const yMin = yOffset + trackH - paddingH - ((sMin - minVal) / valRange) * drawH;
      const yMax = yOffset + trackH - paddingH - ((sMax - minVal) / valRange) * drawH;

      if (px === 0) {
        ctx.moveTo(px, yMin);
        ctx.lineTo(px, yMax);
      } else {
        ctx.lineTo(px, yMin);
        ctx.lineTo(px, yMax);
      }
    }
    ctx.stroke();

    if (isVol) {
      ctx.fillStyle = color;
      this.breaths.forEach(b => {
        if (b.peak_s >= this.viewMinS && b.peak_s <= this.viewMaxS) {
          const px = ((b.peak_s - this.viewMinS) / duration) * w;
          const py = yOffset + trackH - paddingH - ((b.vt - minVal) / valRange) * drawH;
          ctx.beginPath();
          ctx.arc(px, py, 4, 0, Math.PI * 2);
          ctx.fill();
        }
      });
    }

    ctx.fillStyle = 'rgba(255, 255, 255, 0.85)';
    ctx.font = '600 11px Inter';
    ctx.fillText(label, 12, yOffset + 20);
  }

  drawOverlaySignal(ctx, arr, yOffset, trackH, color, alpha, isVol = false) {
    const rect = this.chartViewport.getBoundingClientRect();
    const w = rect.width;
    const range = this.getVisibleIndices();

    let minVal = Infinity, maxVal = -Infinity;
    for (let i = range.start; i <= range.end; i++) {
      if (arr[i] < minVal) minVal = arr[i];
      if (arr[i] > maxVal) maxVal = arr[i];
    }
    const valRange = (maxVal - minVal) || 1;
    const drawH = trackH * 0.8;

    ctx.strokeStyle = color;
    ctx.globalAlpha = alpha;
    ctx.lineWidth = 1.8;
    ctx.beginPath();

    const duration = this.viewMaxS - this.viewMinS;
    const numPixels = Math.ceil(w);
    for (let px = 0; px < numPixels; px++) {
      const t0 = this.viewMinS + (px / numPixels) * duration;
      const t1 = this.viewMinS + ((px + 1) / numPixels) * duration;

      const idx0 = this.findTimeIndex(t0);
      const idx1 = Math.max(idx0, this.findTimeIndex(t1));

      let sMin = arr[idx0], sMax = arr[idx0];
      for (let j = idx0; j <= idx1 && j < arr.length; j++) {
        if (arr[j] < sMin) sMin = arr[j];
        if (arr[j] > sMax) sMax = arr[j];
      }

      const yMin = yOffset + trackH * 0.9 - ((sMin - minVal) / valRange) * drawH;
      const yMax = yOffset + trackH * 0.9 - ((sMax - minVal) / valRange) * drawH;

      if (px === 0) {
        ctx.moveTo(px, yMin);
        ctx.lineTo(px, yMax);
      } else {
        ctx.lineTo(px, yMin);
        ctx.lineTo(px, yMax);
      }
    }
    ctx.stroke();
    ctx.globalAlpha = 1.0;
  }

  findTimeIndex(t) {
    const n = this.time_s.length;
    if (!n) return 0;
    let idx = Math.round((t / (this.time_s[n - 1] || 1)) * (n - 1));
    return Math.max(0, Math.min(n - 1, idx));
  }

  getVisibleIndices() {
    let start = 0, end = this.time_s.length - 1;
    for (let i = 0; i < this.time_s.length; i++) {
      if (this.time_s[i] >= this.viewMinS) {
        start = Math.max(0, i - 1);
        break;
      }
    }
    for (let i = this.time_s.length - 1; i >= 0; i--) {
      if (this.time_s[i] <= this.viewMaxS) {
        end = Math.min(this.time_s.length - 1, i + 1);
        break;
      }
    }
    return { start, end };
  }

  handleMouseMove(e, canvas) {
    if (this.isDragging) {
      const rect = this.chartViewport.getBoundingClientRect();
      const dx = e.clientX - this.dragStartX;
      const dt = -(dx / rect.width) * (this.dragStartMaxS - this.dragStartMinS);
      const span = this.dragStartMaxS - this.dragStartMinS;

      const totalMax = this.time_s[this.time_s.length - 1] || 60.0;
      let newMin = this.dragStartMinS + dt;
      let newMax = this.dragStartMaxS + dt;

      if (newMin < 0) {
        newMin = 0;
        newMax = span;
      }
      if (newMax > totalMax) {
        newMax = totalMax;
        newMin = totalMax - span;
      }

      this.viewMinS = newMin;
      this.viewMaxS = newMax;
      this.renderCharts();
      return;
    }

    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const frac = x / rect.width;
    const targetT = this.viewMinS + frac * (this.viewMaxS - this.viewMinS);

    const bestIdx = this.findTimeIndex(targetT);
    const t = this.time_s[bestIdx].toFixed(2);
    const fVal = this.flow[bestIdx].toFixed(1);
    const cVal = this.ch4[bestIdx].toFixed(1);
    const vVal = this.vol[bestIdx].toFixed(1);

    this.tooltipBar.innerHTML = `
      Idő: <strong>${t} s</strong> &nbsp;|&nbsp;
      1. csatorna (áramlás): <strong style="color:var(--color-flow)">${fVal} ml/s</strong> &nbsp;|&nbsp;
      4. csatorna (CO2): <strong style="color:var(--color-co2)">${cVal} Hgmm</strong> &nbsp;|&nbsp;
      Volumen: <strong style="color:var(--color-vol)">${vVal} ml</strong>
    `;

    this.renderCharts();
    const ctx = canvas.getContext('2d');
    ctx.strokeStyle = 'rgba(0, 240, 255, 0.5)';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, canvas.height);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  handleMouseDown(e) {
    this.isDragging = true;
    this.dragStartX = e.clientX;
    this.dragStartMinS = this.viewMinS;
    this.dragStartMaxS = this.viewMaxS;
  }

  handleMouseUp() {
    this.isDragging = false;
  }

  handleWheel(e) {
    e.preventDefault();
    const zoomFactor = e.deltaY > 0 ? 1.15 : 0.87;
    const span = (this.viewMaxS - this.viewMinS) * zoomFactor;
    const minSpan = 2.0;
    const totalMax = this.time_s[this.time_s.length - 1] || 60.0;
    const maxSpan = totalMax;

    const clampedSpan = Math.max(minSpan, Math.min(maxSpan, span));
    const center = (this.viewMinS + this.viewMaxS) / 2.0;
    let newMin = center - clampedSpan / 2.0;
    let newMax = center + clampedSpan / 2.0;

    if (newMin < 0) {
      newMin = 0;
      newMax = clampedSpan;
    }
    if (newMax > totalMax) {
      newMax = totalMax;
      newMin = totalMax - clampedSpan;
    }

    this.viewMinS = newMin;
    this.viewMaxS = newMax;
    this.renderCharts();
  }

  exportSignalsCSV() {
    if (!this.time_s.length) return;
    let csv = "time_s,respiratory_flow_mls,time_capnogram_mmhg,calculated_volume_ml\n";
    for (let i = 0; i < this.time_s.length; i++) {
      csv += `${this.time_s[i].toFixed(4)},${this.flow[i].toFixed(2)},${this.ch4[i].toFixed(2)},${this.vol[i].toFixed(2)}\n`;
    }
    this.downloadFile(csv, `${this.rawFile.filename}_kalibralt_jelek.csv`, "text/csv");
  }

  exportBreathsCSV() {
    if (!this.breaths.length) return;
    let csv = "cycle_index,start_s,peak_s,end_s,tidal_volume_ml,vd_fowler_ml,peco2_mmhg,paco2_mmhg,bohr_dead_space_pct,bohr_dead_space_ml,vco2_ml,ti_s,te_s,ie_ratio,respiratory_rate_bpm,etco2_mmhg\n";
    this.breaths.forEach(b => {
      csv += `${b.index},${b.start_s.toFixed(3)},${b.peak_s.toFixed(3)},${b.end_s.toFixed(3)},${b.vt.toFixed(2)},${b.fowler_vd.toFixed(1)},${b.peco2.toFixed(2)},${b.paco2.toFixed(2)},${(b.bohr_ratio*100).toFixed(1)},${b.bohr_vd.toFixed(1)},${b.vco2_ml.toFixed(2)},${b.ti.toFixed(3)},${b.te.toFixed(3)},"1:${b.ieRatioVal.toFixed(2)}",${b.rr.toFixed(1)},${b.et_co2.toFixed(1)}\n`;
    });
    this.downloadFile(csv, `${this.rawFile.filename}_volumetrikus_kapnografia.csv`, "text/csv");
  }

  downloadFile(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }
}

window.addEventListener("DOMContentLoaded", () => {
  window.respViewer = new RespViewerApp();
});
