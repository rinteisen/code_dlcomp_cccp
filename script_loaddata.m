clear;
T = 80;
s_set = [1,2,3,4,5,7,8,9,10];
for s = s_set
	data(s,1) = load(['./proposed/result_UE_20_CH_6_round_',int2str(s)]);
	data(s,2) = load(['./fmincon_pa/result_UE_20_CH_6_round_',int2str(s)]);
	data(s,3) = load(['./fmincon_pa_1/result_UE_20_CH_6_round_',int2str(s)]);
end

propfair = zeros(3,T);
time = zeros(3,1);
for m = 1:3
	a = vec2mat([data(:,m).proportional_fair_index],80);
	b = vec2mat([data(:,m).time],80);
	propfair(m,:) = mean(a,1);
	time(m) = mean(b(:,end));
end

figure;
hold on;
plot(propfair(1,:),'r');
plot(propfair(2,:),'g');
plot(propfair(3,:),'b');
hold off;
grid on
legend('(RA,PA) = (proposed,proposed)','(RA,PA) = (proposed,cccp)','(RA,PA) = (proposed,cccp1)');
xlabel('scheduling time slot');
ylabel('proportional fairness');
time

% x = 1:80;

% figure
% grid on
% [AX,H1,H2] = plotyy(x,[propfair(1,:);propfair(2,:);propfair(3,:)],...
% 					x,[exp((propfair(2,:)-propfair(1,:))/20);exp((propfair(3,:)-propfair(1,:))/20)]);
% axes(AX(1));
% ylabel('proportional fair index');
% set(AX(1),'Ylim',[-70,10]);
% set(AX(1),'yTick',[-70:10:10]);
% axes(AX(2));
% set(AX(2),'Ylim',[0.75,1]);
% set(AX(2),'yTick',[0.75:0.05:1])
% ylabel('ratio of average rate per UE');
% legend([H1;H2],'w/ cooperation','w/o cooperation');
% xlabel('scheduling time slot');
% set(AX,'FontSize',14);