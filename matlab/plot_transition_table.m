%% PLOT ALL POSSIBLE REGION TRANSITIONS
figure
%title('Possible region transitions')
xlabel('From region','interpreter','latex')
ylabel('To region','interpreter','latex')
set(gca,'fontsize',24);

hold on
nr_of_regions = 2*m^(n-1);
 for k = 1:nr_of_regions
    for c = 1:(size(Reachable_regions_regDetQ{k,1},2))
        region_to = Reachable_regions_regDetQ{k,1}(c);
        plot(k,region_to,'k*')
    end
     
 end
 drawnow
 