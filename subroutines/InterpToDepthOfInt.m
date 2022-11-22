function [var_int]=InterpToDepthOfInt(var,depth_levels,depth_of_int,ctr)

% Kori-ULB
% Interpolation of ocean variables to ice shelf depth

    var_int=zeros(ctr.imax,ctr.jmax);    
    var_up=var(:,:,1);
    var_int(depth_of_int>=depth_levels(1))=var_up(depth_of_int>=depth_levels(1));
    var_down=var(:,:,end);
    var_int(depth_of_int<=depth_levels(end))=var_down(depth_of_int<=depth_levels(end));
    for k=length(depth_levels):-1:2
        var_up=var(:,:,k);
        var_down=var(:,:,k-1);
        var_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))=...
            ((depth_levels(k)-depth_of_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))).*var_down(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))...
            +(depth_of_int(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1))-depth_levels(k-1)).*var_up(depth_of_int>=depth_levels(k) & depth_of_int<=depth_levels(k-1)))...
            ./(depth_levels(k)-depth_levels(k-1));
    end
%     var_int(icemask_shelves<5.e-1 | depth_of_int>0.)=-9999.9;

end


