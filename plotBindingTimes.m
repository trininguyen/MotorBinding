clear all
format compact

dist = [0 20 40 60];
longDist = [0 20 40 60 80];
longerDist = [0 20 40 60 80 100];
Motlength = [33.1429 45.7143 60.8980];
k_ADP_off = 0.008;
k_ADP_on = 883.8584;
k_ADP_fast_off = 2.12;
k_MT_dir_on = 70.652;
k_MT_dir_off = 0.2;
D_m = [100 1994 10000];
k_weak = 0.002;
car_radi = [280];

%experimental data
expMean = [1.40181 2.44253 6.90707 17.69813 0;...
     1.01732 1.34771 2.61329 5.32338 0;...
     1.01553 0.9931 1.10413 2.1596 10.12818];
expSEM = [0.09242 0.19848 1.00077 2.90438 0;...
     0.07047 0.09428 0.24651 0.67904 0;...
     0.06051 0.07193 0.09975 0.28432 1.88437];

c = hsv(length(D_m)); %or whatever parameter we want to sweep

for k = 1:length(Motlength)
    figure(k);clf;hold on
    if k == 3
        dist = longerDist;
    end
    d = dist;
    for x = 1:length(D_m) %or whatever parameter we want to sweep
        count = 1;
        for i = d
            parfor j = 1:800
				%for Simple Diffusion model, change "DiffwADPState" to "SimpleDiffusion"
                [timetobind(j)] = DiffwADPState(i,k_ADP_off,k_ADP_on,k_ADP_fast_off,...
                    k_MT_dir_on,k_MT_dir_off,k_CTT_on,k_CTT_off,k_CTT_MT,k_MT_CTT,Motlength(k),D_m(x),k_weak,car_radi);
            end
            timetobind(timetobind>=99) = [];
            MeanTime(x,count) = mean(timetobind);
            SemTime(x,count) = std(timetobind)/sqrt(length(timetobind));
            count = count + 1;
        end
        h(x)=errorbar(d,MeanTime(x,:),SemTime(x,:),'.-','color',c(x,:),'Linewidth',2);
%         if k == 3
%             ex=errorbar(dist,expMean(k,:),expSEM(k,:),'k.-','markersize',30,'Linewidth',2);
%         else
%             ex=errorbar(dist,expMean(k,1:end-1),expSEM(k,1:end-1),'k.-','markersize',30,'Linewidth',2);
%         end
    end

    xlabel('Distance between Cargo and MT (nm)','color','k')
    ylabel('Mean Time to Bind (secs)','color','k')
    title(['Motor Length ' num2str(Motlength(k)) ' (nm)'],'color','k');
    box on
    grid on
    set(gca,'fontsize',15)
    leg = legend([h],[string(D_m)],'location','best'); %or whatever parameter we want to sweep
    title(leg,'Cargo Radius (nm)','color','k')
    xlim([-1 dist(end)+1])
end
     
