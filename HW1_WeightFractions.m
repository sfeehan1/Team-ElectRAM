
%% L/D Ratio Calculation
function [LD_ratio] = compute_LD_ratio(q, CD0, W, W_TO, W_S, e, AR)
W_by_W_TO = W/W_TO;
W_by_S = W_by_W_TO*W_S;
LD_ratio = 1 / ((q*CD0/(W_by_S))+((W_by_S)/(q*pi*e*AR)));
end

%% Weight Fraction Calculation
function [WF] = compute_weightfraction(TSFC, R, Vend, LD_ratio)
WF = exp(-((R*TSFC)/(Vend*LD_ratio)));
end

%% Climb segment
function [W_out, fuel_used] = segment_climb(W_TO, W_in, Mach)
WF_Climb = 1.0065 - 0.0325 * Mach;
fuel_used = (1-WF_Climb) * W_in;
W_out = W_in - fuel_used;
end

%% Combat Segment
function [W_out, fuel_used] = segment_combat(W_in, time, TSFC, payload)
fuel_used = time * 9906.98 * TSFC;
W_out = W_in - fuel_used - payload;
end

%% Cruise Segment
function [W_out, fuel_used] = segment_cruise(W_in, W_S, TSFC, Distance, Mach, a, q, CD0, e, AR, W_TO)
V = Mach * a;
LD = compute_LD_ratio(q, CD0, W_in, W_TO, W_S, e, AR);
WF_Cruise = compute_weightfraction(TSFC, Distance, V, LD);
fuel_used = W_in * (1 - WF_Cruise);
W_out = W_in - fuel_used;
end

%% Dash Segment
function [W_out, fuel_used] = segment_dash(W_in, W_TO, W_S, q, CD0, e, AR, TSFC, Distance, V)
LD = compute_LD_ratio(q, CD0, W_in, W_TO, W_S, e, AR);
WF_Dash = compute_weightfraction(TSFC, Distance, V, LD);
fuel_used = W_in * (1 - WF_Dash);
W_out = W_in - fuel_used;
end

%% Landing Segment
function [W_out, fuel_used] = segment_landing(W_in, W_TO)
    WF_Landing = 0.995;
    fuel_used = W_in * (1 - WF_Landing);
    W_out = W_in - fuel_used;
end

%% Lioter Segment
function [W_out, fuel_used] = segment_loiter(W_TO, W_in, W_S, q, CD0, e, AR, time, TSFC)
    LD = compute_LD_ratio(q, CD0, W_in, W_TO, W_S, e, AR);
    WF_Loiter = exp(-(time*60*TSFC/LD));
    fuel_used = W_in * (1 - WF_Loiter);
    W_out = W_in - fuel_used;
end

%% Takeoff Segment
function [W_out, fuel_used] = segment_takeoff(W_in)
WF = 0.95;
W_out = W_in * WF;
fuel_used = W_in - W_out;
end
