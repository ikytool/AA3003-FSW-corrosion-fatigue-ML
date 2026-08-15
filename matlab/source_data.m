function S = source_data()
%SOURCE_DATA Published constants from the two source papers.
%   S = source_data() returns every number taken from:
%     [JMRT]     Chekalil et al., J. Mater. Res. Technol. 33 (2024) 2353-2364
%     [COLLOIDS] Chekalil et al., Colloids Surf. A 680 (2024) 132673
%   Values marked VERIFY were read from the PDF text layer and must be
%   checked against the typeset tables before submission.

% ---------------- Fatigue: published Paris fits [JMRT Figs. 12/14] ----
% da/dN = C * dK^m, da/dN in mm/cycle, dK in MPa*sqrt(m)   (VERIFY units)
S.paris.condition = {'BM', 'FSW', 'FSW.C 25C', 'FSW.C 55C', 'FSW.C 85C'};
S.paris.C  = [3e-11, 9e-12, 7e-12, 2e-11, 7e-12];
S.paris.m  = [3.7741, 4.246, 4.5757, 4.0869, 5.0928];
S.paris.R2 = [0.993, 0.9789, 0.987, 0.9743, 0.9835];
S.paris.dK_range = [3.75 18];           % MPa*sqrt(m)
S.paris.T_C = [NaN, NaN, 25, 55, 85];   % pre-corrosion temperature, degC

% Fatigue test parameters [JMRT]
S.test.R = 0.1; S.test.freq_Hz = 15; S.test.Fmax_kN = 3.5;
S.test.thickness_mm = 1.8; S.test.precorrosion_h = 6;

% Low-dK inhibitor efficiency endpoints reported in [JMRT] (percent)
S.fatigue.eta_T_C    = [25 55 85];
S.fatigue.eta_lowdK  = [8.2 9.6 48.5];   % VERIFY against Fig. 14 text

% Service-life reductions of FSW.C vs blank FSW reported in [JMRT]
S.fatigue.life_drop_55C_pct = 2.42;
S.fatigue.life_drop_85C_pct = 61.64;

% ---------------- Gravimetry [COLLOIDS] ------------------------------
% Corrosion rate W vs temperature, solution Y (0.1 M HCl + 35 g/L NaCl).
% The published unit label (g.mm^-2.h^-1) is a misprint; a dimensional
% audit of the printed mass losses (~0.004 g), specimen area
% (30 x 27 x 1.8 mm, ~1.8e-3 m^2) and 3 h immersion gives g.m^-2.h^-1.
% Only rate RATIOS enter the analysis, so Ea is unaffected.
S.grav.T_C     = [25 45 65 85];
S.grav.W_BM_Y  = [0.78882 2.26503 5.06672 9.07076];
% FSW series in solution Y: only the 85 degC value is unambiguous in the
% text layer (Table 4, 0% inhibitor column). Others: VERIFY from Table 3.
S.grav.W_FSW_Y = [NaN NaN NaN 9.1166];

% Table 4: effect of ethylene glycol concentration at 85 degC, solution Y
S.grav.eg_conc_pct = [0 10 20 30 40];
S.grav.W_BM_eg  = [9.07076 5.14374 4.97412 4.27920 3.40393];
S.grav.W_FSW_eg = [9.1166  4.69528 4.60503 4.24944 2.97671];
S.grav.eta_BM_eg  = [NaN 43.2931 45.1631 52.8242 62.4736];  % percent
S.grav.eta_FSW_eg = [NaN 48.2372 49.2321 53.1523 67.1834];  % percent

% Table 5: activation energies AS PUBLISHED (J/mol as printed; the
% analysis re-derives Ea from the corrosion rates and quotes these for
% traceability).  Order: [BM-Y, FSW-Y, BM-Z, FSW-Z]
S.grav.Ea_published_Jmol = [605.751921 643.000282 685.14373 755.406864];
S.grav.Ea_published_R2   = [0.977 0.986 0.916 0.907];

% Table 7: Langmuir adsorption (solution Z)
S.grav.Kads_BM = 7.56029;  S.grav.dG_BM_kJmol = -14.972;
S.grav.Kads_FSW = 8.87000; S.grav.dG_FSW_kJmol = -15.369;

% ---------------- Constants ------------------------------------------
S.Rgas = 8.314;                          % J/(mol K)
S.TK   = @(TC) TC + 273.15;              % degC -> K
end
