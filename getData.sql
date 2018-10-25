SELECT TOP 250000
  p.modelMag_u,
  p.modelMag_g,
  p.modelMag_r,
  p.modelMag_i,
  p.modelMag_z,
  p.extinction_u,
  p.extinction_g,
  p.extinction_r,
  p.extinction_i,
  p.modelMag_u - 5.*log(299792.458*s.z*1000000./72.5) - p.extinction_u as M_u,
  p.modelMag_g - 5.*log(299792.458*s.z*1000000./72.5) - p.extinction_g as M_g,
  p.modelMag_r - 5*log(299792.458*s.z*1000000./72.5) - p.extinction_r as M_r,
  p.modelMag_i - 5*log(299792.458*s.z*1000000./72.5) - p.extinction_i as M_i,
  p.modelMag_z - 5.*log(299792.458*s.z*1000000./72.5) - p.extinction_z as M_z,
  s.z,
  p.probPSF
FROM
  PhotoObj AS p
  JOIN SpecObj AS s ON s.bestobjid = p.objid
WHERE
  s.z > 0. AND
  p.probPSF BETWEEN 0.85
  AND 1

