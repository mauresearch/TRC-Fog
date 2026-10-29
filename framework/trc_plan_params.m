function [Pp, bg] = trc_plan_params(ctx, model)
% Planning parameters and background load of a method.
Pp = ctx.P;
Pp.modelThrottle = model.throttle;
if ~model.sleepPlan
    Pp.sleepFrac = 1;
end
if strcmp(model.forecast, 'last')
    bg = ctx.bgLast;
else
    bg = ctx.bgARSENAL;
end
end
