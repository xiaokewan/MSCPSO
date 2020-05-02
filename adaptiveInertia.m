function y = adaptiveInertia(x,favg,wmax,wmin,fmin)	
	%%
    if x >= favg
    %stage1
    w = wmax;

    elseif favg/2<=x<favg
    %stage2
    w = wmin+(x-fmin)*(wmax-wmin)/((favg/2)-fmin);

    elseif x <favg/2
    %stage3
    w = wmax-it*((wmax-wmin)/MaxIt);
    end
    
    y=w;%%
end
