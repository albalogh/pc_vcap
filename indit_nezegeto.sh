#!/usr/bin/env bash
# CapnoView indítóparancsfájl
# Elindítja a helyi Python szervert (8088-as port), és megnyitja a böngészőt az interaktív nézegetővel.

PORT=8088
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/resp_viewer" && pwd)"

# Keressük meg a megfelelő python3 futtatót
PYTHON_CMD="python3"
if [ -x "/home/balogh/Dokumentumok/lc/.venv/bin/python3" ]; then
    PYTHON_CMD="/home/balogh/Dokumentumok/lc/.venv/bin/python3"
fi

echo "============================================================="
echo "  CapnoView"
echo "  Szerver indítása a(z) ${PORT}-as porton..."
echo "  Fájlok keresési mappája: $(dirname "$APP_DIR")"
echo "============================================================="

# Régi folyamat leállítása, ha fut a 8088-as porton, hogy mindig a legfrissebb kód futhasson
if ss -tulpn 2>/dev/null | grep -q ":${PORT} "; then
    echo "A korábbi szerver leállítása a ${PORT}-as porton..."
    fuser -k -n tcp "$PORT" 2>/dev/null
    sleep 1
fi

# Böngésző megnyitása a háttérben 1.5 másodperc múlva
(
    sleep 1.5
    xdg-open "http://localhost:${PORT}" 2>/dev/null &
) &

# Szerver indítása (setsid a háttérben vagy előtérben futtatáshoz)
cd "$APP_DIR" || exit 1
exec "$PYTHON_CMD" -u app.py -p "$PORT"
