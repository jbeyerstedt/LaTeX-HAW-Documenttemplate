%----------------------------------------------------------------
% dfp_wpp_mes.m
%
%                   DAVIDSON-FLETCHER-POWELL-VERFAHREN
%                     approximiert eine Zielfunktion

function xkplus1=dfp_wpp_mes(  freq, yWunsch, Hfunkt_mfile, Hgradfunkt_mfile, x0, Lp, error_scaling);

% Approximationsparameter festlegen
eps_GoldenSection = 1.e-5;      % eps Goldener Schnitt
max_iter = 10*length(x0);       % MAX_Anzahl_Iterationen
stopbed = 0;                    % Abbruchbedingung, 0 heißt: mache MAX_Anzahl_Iterationen

%---------------- WHILE loop ----------------------------
lpcnt = 0;
while lpcnt <= max_iter
	lpcnt = lpcnt+1;
    f_omc = feval(Hfunkt_mfile, xk, freq, numlen);
    F_c = (sum((f_omc - yWunsch).^Lp))^(1/Lp);
    gk = ((f_omc.' - yWunsch.') ./ F_c).^(Lp-1) * grad;
    gk = gk.';
    sk = Qu * gk;

    it = 0;
    while it <= max_iter
        alphai = alf_Schrittweite_1D * 2^(it-1);

        f_omc = feval(Hfunkt_mfile, (xk + alphai * sk), freq, numlen);
        h_alphai = (sum((f_omc - yWunsch).^Lp))^(1/Lp);

        alpha_iminus2 = alpha_iminus1;
        alpha_iminus1 = alphai;

        if h_alphai > h_alpha_iminus1
            break
        end

        h_alpha_iminus1 = h_alphai;
        it = it+1;
    end

% Fletcher-Powell, Entenmann S. 74, Skript Formel (1)
    Qu = Qu - (Qu * qk * qk.' * Qu) / (qk.' * Qu * qk) + (pk * pk.') / (qk.' * qk);
% Fletcher-Broyden, Entenmann S. 79
%  	Qu=Qu-...

% mit neuen Werten nächste Iteration berechnen....
    xk = xkplus1;

end; % WHILE lpcnt=1:r...
