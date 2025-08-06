function [flag]=NaNcheck(flag,H,ux,uy)

% Kori-ULB
% Check at any time step whether NaN or unusual values are occuring to stop
% the run before crashing

    if any(isnan(H(:)))
        fprintf('\n\n NaNs occured in H\n');
        flag=true;
    end
    if any(isnan(ux(:)))
        fprintf('\n\n NaNs occured in ux\n');
        flag=true;
    end
    if any(isnan(uy(:)))
        fprintf('\n\n NaNs occured in uy\n');
        flag=true;
    end
    if any(isinf(H(:)))
        fprintf('\n\n Inf occured in H\n');
        flag=true;
    end
    if any(isinf(ux(:)))
        fprintf('\n\n Inf occured in ux\n');
        flag=true;
    end
    if any(isinf(uy(:)))
        fprintf('\n\n Inf occured in uy\n');
        flag=true;
    end
    if any(H(:)>1e4)
        fprintf('\n\n A unrealistic large ice thickness (H>1e4) occured\n');
        flag=true;
    end

end


