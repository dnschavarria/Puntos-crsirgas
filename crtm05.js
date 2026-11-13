// Conversión geográficas (lat/lon) <-> CR-SIRGAS / CRTM05 (EPSG:8908)
// Transversa de Mercator, serie de Krüger (orden 6, error < 1 mm en Costa Rica).
// Parámetros: lat0 = 0°, lon0 = -84°, k0 = 0.9999, FE = 500 000 m, FN = 0 m, elipsoide GRS80.
// Se asume WGS84 ≈ CR-SIRGAS (transformación nula, como EPSG:9752), válido para
// la precisión de un GNSS de celular (metros).
(function (root) {
  const a = 6378137.0;
  const f = 1 / 298.257222101; // GRS80
  const k0 = 0.9999, lon0 = -84 * Math.PI / 180, FE = 500000, FN = 0;
  const n = f / (2 - f);
  const n2 = n * n, n3 = n2 * n, n4 = n3 * n, n5 = n4 * n, n6 = n5 * n;
  const A = a / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256);
  const alpha = [0,
    n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
    13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
    61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
    49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
    34729 * n5 / 80640 - 3418889 * n6 / 1995840,
    212378941 * n6 / 319334400];
  const beta = [0,
    n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
    n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
    17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
    4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
    4583 * n5 / 161280 - 108847 * n6 / 3991680,
    20648693 * n6 / 638668800];
  const e = Math.sqrt(f * (2 - f));

  function toCRTM05(latDeg, lonDeg) {
    const phi = latDeg * Math.PI / 180, lam = lonDeg * Math.PI / 180 - lon0;
    const t = Math.sinh(Math.atanh(Math.sin(phi)) - e * Math.atanh(e * Math.sin(phi)));
    const xi1 = Math.atan2(t, Math.cos(lam));
    const eta1 = Math.atanh(Math.sin(lam) / Math.sqrt(1 + t * t));
    let xi = xi1, eta = eta1;
    for (let j = 1; j <= 6; j++) {
      xi += alpha[j] * Math.sin(2 * j * xi1) * Math.cosh(2 * j * eta1);
      eta += alpha[j] * Math.cos(2 * j * xi1) * Math.sinh(2 * j * eta1);
    }
    return { E: FE + k0 * A * eta, N: FN + k0 * A * xi };
  }

  function fromCRTM05(E, N) {
    const xi = (N - FN) / (k0 * A), eta = (E - FE) / (k0 * A);
    let xi1 = xi, eta1 = eta;
    for (let j = 1; j <= 6; j++) {
      xi1 -= beta[j] * Math.sin(2 * j * xi) * Math.cosh(2 * j * eta);
      eta1 -= beta[j] * Math.cos(2 * j * xi) * Math.sinh(2 * j * eta);
    }
    const chi = Math.asin(Math.sin(xi1) / Math.cosh(eta1));
    const lam = Math.atan2(Math.sinh(eta1), Math.cos(xi1));
    // latitud a partir de la conforme (iteración)
    let phi = chi;
    for (let i = 0; i < 10; i++) {
      const s = Math.sin(phi);
      phi = 2 * Math.atan(Math.tan(Math.PI / 4 + chi / 2) *
        Math.pow((1 + e * s) / (1 - e * s), e / 2)) - Math.PI / 2;
    }
    return { lat: phi * 180 / Math.PI, lon: (lam + lon0) * 180 / Math.PI };
  }

  const api = { toCRTM05, fromCRTM05, EPSG: 8908,
    PROJ: '+proj=tmerc +lat_0=0 +lon_0=-84 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs' };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.CRTM05 = api;
})(this);
