%% Plot Figures
%  This script plots the results obtained from the simulation

%% Plot ETC triggering sequence
figure
plot(trig_time,'*','MarkerSize',6)
%axis([0 520 0 0.35])
%set(gca,'fontsize',30,'fontname','Arial');
set(gca,'fontsize',24);
%title('Triggering time sequence of the system')
xlabel('Trigger instance','interpreter','latex')
ylabel('Inter-sample time (s)','interpreter','latex')
hold on,
plot((1:size(trig_time,2)),Min_reg_seq,'k-',(1:size(trig_time,2)),Max_reg_seq,'r-','LineWidth',1.4)
%legend('ETC triggering','$\underline{\tau}_s$','$\overline{\tau}_s$','interpreter','latex')
 
[~,hObj] = legend('ETC triggering','$\underline{\tau}_s$','$\overline{\tau}_s$');                   % return the handles array
 hL=findobj(hObj,'type','line');  % get the lines, not text
 hT=findobj(hObj,'type','text');  % get the text
 set(hL,'linewidth',1)            % set their width property
 set(hT,'interpreter','latex','fontsize',30);

%% Plot measurement error evolution

% figure
% hold on
% plot((0:ts:time_end-ts),alpha^(0.5)*x_norm_resp,'LineWidth',1.4)
% set(gca,'fontsize',30,'fontname','Arial');
% plot((0:ts:time_end-ts),e_norm_resp,'r','LineWidth',1.4)
% title('Evolution of $\alpha \times |x|$ and $|e|$','interpreter','latex')
% xlabel('Time ($s$)','interpreter','latex')
% ylabel('$\alpha |x|$ and $|e|$','interpreter','latex')
% hleg = legend('$\alpha |x|$','$|e|$');
% set(hleg,'interpreter','latex');
%% Plot state evolution during simulation
% NOTE: adjust legend depending on which control system is simulated!!!
% SIMULATION
figure
hold on
plot((0:ts:time_end-ts),x_resp(1,:),'b')
plot((0:ts:time_end-ts),x_resp(2,:),'r')
% plot((0:ts:time_end-ts),x_resp(3,:),'k')
%plot((0:ts:time_end-ts),x_resp(4,:),'g')

% if n > 4
%     plot((0:ts:time_end-ts),x_resp(4,:),'g')
% end

set(gca,'fontsize',24,'fontname','Arial');
%title('Evolution of states')
xlabel('Time $(s)$','interpreter','latex')
ylabel('$x_i$','interpreter','latex')
%axis([0 time_end -6 6.5])
%legend('Distance error $(m)$','Velocity error $(m/s)$','Vehicle acceleration $(m/s^2)$','interpreter','latex')

[~,hObj] = legend('$x_1-x_2 (m)$','$\dot{x_1}(m/s)$','$x_2 - z_r (m)$','$\dot{x_2}(m/s)$'); % return the handles array
hL=findobj(hObj,'type','line');  % get the lines, not text
hT=findobj(hObj,'type','text');  % get the text
set(hL,'linewidth',2)            % set their width property
set(hT,'interpreter','latex','fontsize',30);

%% Plot measurement errors

% figure
% hold on
% plot((0:ts:time_end-ts),e_resp(1,:),'b')
% plot((0:ts:time_end-ts),e_resp(2,:),'r')
% plot((0:ts:time_end-ts),e_resp(3,:),'g')
% title('State Errors')
% xlabel('time ($s$)','interpreter','latex')
% ylabel('$e_i$','interpreter','latex')
% hleg = legend('$e_1$','$e_2$','$e_3$','$e_4$');
% set(hleg,'interpreter','latex');

%% Find maximum difference between upper and lower time bound for all regions

max_diff_tau = max(Tau_s_max-Tau_s_opt)

%% Check for impossible transitions (transitions not predicted to be possible by reachability analysis

impossible_transitions = 0;
%trouble_poly = [];             
% To save polyhedra for which reachability analysis failed; requires the
% previous section to be uncommented!!!

for i = 1:(size(reg_Seq,2)-1)
    region = reg_Seq(i);
    region_next = reg_Seq(i+1);
    
    % If a region is not among the reachable regions for the previous
    % region, the transition was not predicted to be possible by the
    % reachability analysis
    if ismember(region_next, Reachable_regions_regDetQ{region,1}) == 0
        impossible_transitions = impossible_transitions+1;
        %trouble_poly = [trouble_poly, reg_Seq_Polyhedra(i)];   
    end
end

% Polyhedron for which reachability analysis is not correct
%trouble_poly = unique(trouble_poly);

% Number of impossible transitions (should be zero)
impossible_transitions

%% PLOT ALL POSSIBLE REGION TRANSITIONS
figure
%title('Possible region transitions')
xlabel('From region','interpreter','latex')
ylabel('To region','interpreter','latex')
set(gca,'fontsize',24);

hold on
 for k = 1:q
     
    for c = 1:(size(Reachable_regions_regDetQ{k,1},2))
        region_to = Reachable_regions_regDetQ{k,1}(c);
        plot(k,region_to,'k*')
    end
     
 end
 