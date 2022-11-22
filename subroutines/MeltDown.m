function MeltDown()

% Kori-ULB
% No Comment

    figure;
    for i=1:20
        k=rand-0.5;
        hold on; clf;
        text(1+5*k,6+5*k,'All is lost','FontSize',24,'color',rand(1,3));
        text(1+5*k,4+5*k,'Complete meltdown','FontSize',24,'color',rand(1,3));
        axis([0 10 0 10]);
        axis off;
        hold off;
        pause(1);
    end
    
end


