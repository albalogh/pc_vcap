# CapnoView

**CapnoView** egy nagy teljesítményű, kutatási és klinikai célra fejlesztett interaktív vizualizációs webalkalmazás és dekódoló rendszer, amely a mainstream időkapnográffal és légzési pneumotachográffal rögzített **`.inp`** (bináris) és **`.csv`** adatsorok elemzésére szolgál.

---

## 🌟 Fő funkciók és jellemzők

1. **`.inp` bináris fájl dekódolás & CSV konverzió (`decode_inp.py` & HTTP API)**
   - 44 bájtos bináris fejléc, méréshatár és Pascal-string metaadatok teljes körű olvasása.
   - 1. csatorna (spontán légzési áramlás) és 4. csatorna (mainstream időkapnogram) közvetlen feldolgozása.

2. **Volumetrikus kapnográfia ($CO_2$ vs. kilégzett volumen - 4. és 5. lépés)**
   - **Fowler-féle anatómiai holttér ($V_{D,Fowler}$):** A volumetrikus kapnogram II. fázisának inflexiós pontjából.
   - **Kilégzett átlagos $CO_2$ ($PECO_2$) & Alveoláris átlagos $CO_2$ ($PACO_2$):** Volumen-súlyozott integrálással és az alveoláris plató fázisának elemzésével.
   - **Bohr-féle élettani holttér ($V_{DB}$):** A Bohr-egyenlet alapján mind százalékos arányban (%), mind abszolút térfogatban (ml):
     $$V_{DB} \text{ (\%)} = \frac{PACO_2 - PECO_2}{PACO_2} \cdot 100$$
   - **Kilégzett $CO_2$-térfogat ($VCO_2$):** Egy ciklus alatt eliminált szén-dioxid mennyiség (ml/ciklus és ml/min).

3. **14 colos kijelzőre optimalizált 50/50 osztott nézet (Tekerés nélkül)**
   - A bal oldali 50%-on található az interaktív grafikonos nézegető (*Egyszerre 3 sávos*, *Kombinált overlay*, *Légzési hurok*, *📈 Volumetrikus kapnogram*).
   - A jobb oldali 50%-on fut a légzésenkénti 14 oszlopos interaktív táblázat; a sorra kattintva a grafikon ráközelít a kiválasztott légzési ciklusra.
   - A képernyő alsó teljes szélességű sávjában (6. lépés) olvasható az összesített élettani statisztika és végezhető el a CSV adatexport.

4. **Utólagos kalibráció (2. lépés)**
   - Közvetlenül paraméterezhető CO2 és áramlás kalibrációs szorzók (`gain`) és alapvonal eltolások (`zero`).

---

## 🚀 Indítás

```bash
bash indit_nezegeto.sh
```

A szerver a **8088**-as porton indul, a felület a böngészőben a [http://localhost:8088](http://localhost:8088) címen érhető el.

---

## 📁 Projekt struktúra

```
pc_vcap/
├── resp_viewer/
│   ├── app.py         # Python backend szerver és volumetrikus kapnográfia számítás
│   ├── app.js         # Frontend logika (50/50 osztott nézet, canvas rendering)
│   ├── index.html     # Felület (1-6. lépés logikai struktúra)
│   └── style.css      # Dark Medical Glassmorphism téma
├── decode_inp.py      # Önálló parancssori dekódoló script (.inp -> .csv)
├── indit_nezegeto.sh  # Szerver indítóparancsfájl
└── README.md
```
