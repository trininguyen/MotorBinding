function [timetobind,MarchovState,tFinal] = DiffwADPState(dist,k_ADP_off,k_ADP_on,k_ADP_fast_off,...
    k_MT_dir_on,k_MT_dir_off,k_CTT_on,k_CTT_off,k_CTT_MT,k_MT_CTT,length,D_m,k_weak,car_radi)

% timetobind = zeros(numSim-1,1);
% j = 1;
% while j < numSim
movie = 0; %if want to visualize, set to 1, else 0
if movie == 1
    capFreq = 1;
    v = VideoWriter('MotorADP.avi');
    v.FrameRate = 10;
    open(v);
end
step_max = 500000; %max number of time steps
step = 1; %time step
time_max = 100; %max number of seconds to simulate
dtmax = 0.004;%0.00305;%0.007; %difference of time between time steps (s)
t = zeros(step_max,1); %time

%car_radi = 280; %cargo radius (nm)
eta = 3e-5;%1.0016e-6;%(pN/nm^2)*s
k_t = 0.19; %stiffness of c terminal tail (pN/nm)
k_m = 0.320; %stiffness of kinesin (pN/nm)
k_b = 0.1895; %stiffness of spring that connects CTT and motor (pN/nm)
%k_s = 100; %stiffness of cargo (pN/nm)
k_laser = [0.045 0.03 0.045];%0.065; %stiffness of laser (pN/nm)
D_t = 2900; %diffusion constant of c terminal tail
%D_m = 2800; %diffusion constant of motor head (nm^2/s)
D_c = 4.114/(6*pi*eta*car_radi); %diffusion constant of cargo 
gamma_m=4.114/D_m;%motor drag coefficient
gamma_t=4.114/D_t;%ctt drag coefficient

l_t = 10.895; %length of C terminal tail (nm)
l_m = length+10; %length of kinesin motor (nm), +10 for binding antibody between motor and cargo (33.1429nm = 406aa, 45.7143nm=560aa, 60.8980=746aa)
numCTT = 130; %number of CTTs in x direction
numP = 6; %number of MT protofilaments to simulate CTTs on, 13 for full circle, 6 for half
CTT_xStart = -500; %x-coordinate of first CTT

%motorTailSpace = 10.8950;%-10; %space between motor head and c terminal tail (nm)

bindRadiusMT = 12.5+2; %radius at which motor head can bind to MT (nm)
bindRadiusCTT = 5; %radius at which motor head can bind to CTT (nm)

% k_ADP_off = 1;
% k_ADP_on = 900;%1262.8; %needs to be less than MT_on
% k_ADP_fast_off = 5;%82.022;
% k_MT_dir_on = 5;%4.3846;
% k_MT_dir_off = 2;
% k_CTT_on = 1000;%39.4616;
% k_CTT_off = 0.01;
% k_CTT_MT = 1500;%31.5692;
% k_MT_CTT = 0.5;

a = [0 12.5+dist 0]; %coordinate of motor anchor
%x_m = [0 l_t+12.5 + motorTailSpace 0]; %coordinate of motor head
x_m_init = [0 a(2)-l_m 0];
if x_m_init(2)<12.5
    x_m_init(2) = 12.5;
end

%a = [0 x_m(2)+l_m 0]; %coordinate of motor anchor
tempAngle = zeros(step_max,1);
theta = a;
c = [0 a(2)+car_radi 0]; %coordinate of cargo center
LenMot=sqrt((a(1)-x_m_init(1))^2+(a(2)-x_m_init(2))^2+(a(3)-x_m_init(3))^2);%length of motor
tempLenMot=LenMot;
closeTails=[];
x_t = zeros(numP+1,3,numCTT);
b = zeros(numP+1,3,numCTT);
distMotorCTT = zeros(numCTT, numP+3);
% for xt = 0:numCTT-1 
%     for yt=-3:numP-3 
%         b(yt+4,:,xt+1)=[(CTT_xStart)+(8*xt) 12.5*cos(12.0830*yt) 12.5*sin(12.0830*yt)];
%         x_t(yt+4,:,xt+1)=[CTT_xStart+(8*xt) b(yt+4,2,xt+1)+l_t*cos(12.0830*yt) b(yt+4,3,xt+1)+l_t*sin(12.0830*yt)];
%         distMotorCTT(xt+1,yt+4) = sqrt((x_m_init(1,1)-x_t(yt+4,1,xt+1))^2 + (x_m_init(1,2)-x_t(yt+4,2,xt+1))^2 + (x_m_init(1,3)-x_t(yt+4,3,xt+1))^2);
%         if distMotorCTT(xt+1,yt+4)<bindRadiusCTT
%             closeTails=[closeTails; xt+1 yt+4];
%         end
%     end
% end


distMotorCTTinit = distMotorCTT;
closeTailsinit = closeTails;

if movie == 1
    figure(1);clf
    r = 12.5;
    [Xc, Yc, Zc] = cylinder(r);
    hc = CTT_xStart*2;Zc = Zc*hc; Zc = Zc - (0.5*hc);
    [X,Y,Z]=sphere(nexttile,100);
    X = X*car_radi; Y=Y*car_radi; Z=Z*car_radi;
    
end
trial = 0;
while step<step_max && t(step)<time_max
%     if trial == 70649
%         movie = 1;
%         capFreq = 1;
%         v = VideoWriter('MotorCTT.avi');
%         v.FrameRate = 10;
%         open(v);
%     end
    if step == 1
        tempTB = 0;
        tempTU = 0;
        timeTailbound = 0;
        trial = trial + 1;
        MotorADPState = 1; %if Motor is in ADP state, set to 1, if in ATP state, set to 2
        tailbound = 0; %if CTT is bound to motor, set to 1
        a = [0 12.5+dist 0]; %coordinate of motor anchor
        c = [0 a(2)+car_radi 0]; %coordinate of cargo center
        MTbound = 0; %if motor is bound to MT, set to 1
        f_mot_head = [0 0 0];
        x_m = zeros(step_max,3);
        x_m(1,:) = x_m_init;
        LenMot=sqrt((a(1)-x_m_init(1))^2+(a(2)-x_m_init(2))^2+(a(3)-x_m_init(3))^2);%length of motor
        tempLenMot=LenMot;
        distMotorCTT = distMotorCTTinit;
        closeTails = closeTailsinit;
        MarchovState = [];
    end
    %Gillespie   
    if MTbound == 0
       if tailbound == 0 && any(any(distMotorCTT <= bindRadiusCTT))
           if MotorADPState == 1
               ADP_off = k_ADP_off;
               ADP_on = 0;
               CTT_on = size(closeTails,1)*k_CTT_on;
               CTT_off = 0;
               MT_on = k_MT_dir_on;
               MT_off = 0;
               MarchovState = [MarchovState 1];
           else
               ADP_off = 0;
               ADP_on = k_ADP_on;
               CTT_on = size(closeTails,1)*k_CTT_on;
               CTT_off = 0;
               MT_on = k_MT_dir_on;
               MT_off = 0;
               MarchovState = [MarchovState 2];
           end
       elseif tailbound==1 && distMotorCTT(boundTailInd(1),boundTailInd(2)) > bindRadiusCTT
           if MotorADPState == 1
               ADP_off = k_ADP_off;
               ADP_on = 0;
               CTT_on = 0;
               CTT_off = k_CTT_off*distMotorCTT(boundTailInd(1),boundTailInd(2));
               MT_on = k_CTT_MT;
               MT_off = 0;
               MarchovState = [MarchovState 3];
           else
               ADP_off = 0;
               ADP_on = k_ADP_on;
               CTT_on = 0;
               CTT_off = k_CTT_off*distMotorCTT(boundTailInd(1),boundTailInd(2));
               MT_on = k_CTT_MT;
               MT_off = 0;
               MarchovState = [MarchovState 4];
           end
       elseif tailbound==1 && distMotorCTT(boundTailInd(1),boundTailInd(2))<=bindRadiusCTT
           CTT_on = 0;
           CTT_off = 0;
           MT_on = k_CTT_MT;
           MT_off = 0;
           if MotorADPState == 1
               ADP_off = k_ADP_off;
               ADP_on = 0;
               MarchovState = [MarchovState 3];
           else
               ADP_off = 0;
               ADP_on = k_ADP_on;
               MarchovState = [MarchovState 4];
           end
       elseif tailbound==0 && all(all(distMotorCTT>bindRadiusCTT))
           if MotorADPState == 1
               ADP_off=k_ADP_off;
               ADP_on = 0;
               CTT_on=0;
               CTT_off=0;
               MT_on=k_MT_dir_on;
               MT_off=0;
               MarchovState = [MarchovState 1];
           else
               ADP_off=0;
               ADP_on = k_ADP_on;
               CTT_on=0;
               CTT_off=0;
               MT_on=k_MT_dir_on;
               MT_off=0;
               MarchovState = [MarchovState 2];
           end
       end
       if sqrt(x_m(step,2)^2+x_m(step,3)^2)>bindRadiusMT
           MT_on = 0;
       end
    else
        if tailbound == 0 && any(any(distMotorCTT <= bindRadiusCTT))
            ADP_off = k_ADP_fast_off;
            ADP_on = 0;
            CTT_on = size(closeTails,1)*k_CTT_on;
            CTT_off = 0;
            MT_on = k_ADP_fast_off;
            MT_off = k_MT_dir_off * exp(norm(f_weak)/4);
            MarchovState = [MarchovState 5];
        elseif tailbound==1 && distMotorCTT(boundTailInd(1),boundTailInd(2)) > bindRadiusCTT
            ADP_off = k_ADP_fast_off;
            ADP_on = 0;
            CTT_on = 0;
            CTT_off = k_CTT_off*distMotorCTT(boundTailInd(1),boundTailInd(2));
            MT_on = k_ADP_fast_off;
            MT_off = k_MT_CTT * sqrt(x_m(step,2)^2+x_m(step,3)^2);
            MarchovState = [MarchovState 6];
        elseif tailbound==1 && distMotorCTT(boundTailInd(1),boundTailInd(2))<=bindRadiusCTT
            ADP_off=k_ADP_fast_off;
            ADP_on=0;
            CTT_on = 0;
            CTT_off=0;
            MT_on=k_ADP_fast_off;
            MT_off=k_MT_CTT * sqrt(x_m(step,2)^2+x_m(step,3)^2);
            MarchovState=[MarchovState 6];
        elseif tailbound==0 && all(all(distMotorCTT>bindRadiusCTT))
            ADP_off=k_ADP_fast_off;
            ADP_on = 0;
            CTT_on = 0;
            CTT_off = 0;
            MT_on = k_ADP_fast_off;
            MT_off = k_MT_dir_off * exp(norm(f_weak)/4);
            MarchovState = [MarchovState 5];
        end
        if MotorADPState==2
            ADP_off = 0;
        end
    end
    array_rates = [ADP_off ADP_on CTT_on CTT_off MT_on MT_off];
    sum_rates = sum(array_rates); %to compute reaction rates
        
    tau = (1/sum_rates)*log(1/rand(1,1)); %how much time before next reaction occurs
    if tau < dtmax
        dt = tau;
        i=1; mu = 0; amu=0; r=rand;
        while amu<r*sum_rates
            mu = mu+1;
            amu = amu + array_rates(i);
            i = i+1;
        end
        
        switch mu
            case 1
                MotorADPState = 2;
            case 2
                MotorADPState = 1;
            case 3
                tailbound = 1;
                tempTB = tempTB + 1;
                timeTBstart = t(step);
                boundTailInd = randi(size(closeTails,1));
                boundTailInd = closeTails(boundTailInd,:);
                %disp('CTT bound');
            case 4
                tailbound = 0;
                tempTU=tempTU + 1;
                timeTailbound = timeTailbound + t(step)-timeTBstart;
            case 5
                MTbound = 1;
                if MotorADPState == 2 %|| finalBind == 1
                    if movie == 1
                        clf;
                        mt = surf(Xc,Yc,Zc);hold on;
                        rotate(mt,[0 1 0], 90)
%                         for xt = 1:numCTT
%                             for yt=1:numP+1
%                                 plot3([x_t(yt,1,xt) b(yt,1,xt)], [x_t(yt,2,xt) b(yt,2,xt)],[x_t(yt,3,xt) b(yt,3,xt)],'color',[0.8500 0.3250 0.0980],'linewidth',2);
%                             end
%                         end
                        plot3(a(step,1),a(step,2),a(step,3),'b.','MarkerSize',20)
                        plot3([x_m(step,1) a(step,1)], [x_m(step,2) a(step,2)],[x_m(step,3) a(step,3)],'linewidth',2);
                        if tailbound==1
                        plot3([x_t(boundTailInd(2),1,boundTailInd(1)) x_m(step,1)],[x_t(boundTailInd(2),2,boundTailInd(1)) x_m(step,2)],[x_t(boundTailInd(2),3,boundTailInd(1)) x_m(step,3)],...
                            'color',[0.3010 0.7450 0.9330],'linewidth',3)
                        end
                        %plot(Cargo_x, Cargo_y);
                        surf(X+c(step,1),Y+c(step,2),Z+c(step,3),'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none')
                        view(2);axis equal
                        xlim([-200 200]);ylim([-12.5 300]);
                        title(['t = ' num2str(t(step)) 's'])
                        drawnow;
                        frame = getframe(gcf);
                        writeVideo(v,frame);
                    end
                    %trial
                    MarchovState = [MarchovState 7];
                    step = step + 1;
                    t(step) = t(step-1) + dt;
                    break;
                else
                   % weak_anchor = [x_m(step,1) x_m(step,2) sqrt((12.5^2)-(x_m(step,2)^2))];
                   P = [x_m(step,2)-0,x_m(step,3)-0];
                    Q = (12.5/norm(P))*P;
                    weak_anchor = [x_m(step,1), Q(1) + x_m(step,2), Q(2) + x_m(step,3)];
                end
            case 6
                MTbound = 0;
                MotorADPState = 1;
                %finalBind = 0;
        end
    else
       dt = dtmax;
    end
    
    %sterics of cargo
%     if sqrt(((-x_m(step,1))^2)+((a(2)+500-x_m(step,2))^2)+((-x_m(step,3)^2)))<=500
% %         f_s = -k_s*(500-norm([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]))*([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]/...
% %             norm([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]));
%         f_s = -k_s*(norm([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]))*([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]/...
%             norm([-x_m(step,1) a(2)+500-x_m(step,2) -x_m(step,3)]));
%     else
%         f_s = 0;
%     end
    
    if tailbound == 1
        if x_m(step,2)<12.5+l_t+bindRadiusCTT && abs(x_m(step,3))<12.5+l_t+bindRadiusCTT
            if round((x_m(step,1)-CTT_xStart)/8) > 2
                num_it_min = round((x_m(step,1)-CTT_xStart)/8)-2;
            else
                num_it_min = 1;
            end
            
            if round((x_m(step,1)-CTT_xStart)/8)<=numCTT-2
                num_it_max = round((x_m(step,1)-CTT_xStart)/8)+2;
            elseif round((x_m(step,1)-CTT_xStart)/8) == numCTT
                num_it_max = round((x_m(step,1)-CTT_xStart)/8)+1;
            else
                num_it_max = round((x_m(step,1)-CTT_xStart)/8);
            end
            
        else
            num_it_min = randi(numCTT);
            num_it_max = num_it_min;
        end
        
        for wt = num_it_min:num_it_max
            for xt=1:numP+1
                if all([wt,xt] == boundTailInd)
                    Px_t_part = x_t(xt,:,wt) + (dt*(-(k_t/gamma_t)*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt)))+k_b*(x_m(step,:)-x_t(xt,:,wt))));
                else
                    Px_t_part = x_t(xt,:,wt) + (-(k_t/gamma_t)*dt*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt))));
                end
                acceptMove = 0;
                ct = 1;
                while acceptMove == 0
                    Px_t = Px_t_part+((sqrt(2*D_t))*sqrt(dt)*randn(1,3)); %Wiener process for c terminal tail
                    %                     if all([wt,xt] == boundTailInd)
                    %                         Px_t = x_t(xt,:,wt) + (dt*(-(k_t/gamma_t)*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt)))+k_b*(x_m(step,:)-x_t(xt,:,wt)))) + ((sqrt(2*D_t))*W_t);
                    %                     else
                    %                         Px_t = x_t(xt,:,wt) + (-(k_t/gamma_t)*dt*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt)))) + ((sqrt(2*D_t))*W_t);
                    %                     end
                    if sqrt(Px_t(1,2)^2+Px_t(1,3)^2) >=12.5 && Px_t(1,2) > -3
                        acceptMove = 1;
                        x_t(xt,:,wt)=Px_t;
                    else
                        ct = ct + 1;
                        if ct >1000
                            break
                        end
                    end
                end
                if ct > 1000
                    disp('Some CTT in CTT-motor bound state cannot move')
                    step = 1;
                    t = zeros(step_max,1);
                    continue
                end
            end
        end
        f_mot_head=-k_m*(LenMot-l_m)*(LenMot>l_m)*((x_m(step,:)-a(step,:))/norm(x_m(step,:)-a(step,:)));
        f_tail_on_head=k_b*(x_t(boundTailInd(2),:,boundTailInd(1))-x_m(step,:));
        c_part = (k_laser.*(c(1,:)-c(step,:))*dt)+((f_mot_head+f_tail_on_head)*dt);
        accept = 0;
        ct = 1;
        while accept == 0
            c(step+1,:)=c(step,:)+(sqrt(dt*2*D_c)*randn(1,3))+c_part;%(k_laser.*(c(1,:)-c(step,:))*dt)+(f_mot_head*dt);
            if c(step+1,2)>12.5+car_radi %&& abs(c(step+1,1))<abs(CTT_xStart)-car_radi-l_m-100
                accept = 1;
            else
                ct = ct+1;
                if ct > 1000
                    break
                end
            end
        end
        if ct>1000
            disp('Cargo cannot move')
            step = 1;
            t = zeros(step_max,1);
            continue
        end
        torque_m = cross((a(step,1)-c(step+1,:)),f_mot_head);
        theta(step+1,:) = theta(step,:) + ((1/(8*pi*eta*(car_radi^3)))*torque_m*dt)+(sqrt(2*(4.114/(8*pi*eta*(car_radi^3)))*dt)*randn(1,3));
        dtheta=theta(step+1,:)-theta(step,:);
        angle = norm(dtheta);
        tempAngle(step) = angle;
        dtheta = dtheta/angle;
        R = [cos(angle)+((dtheta(1)^2)*(1-cos(angle))) (dtheta(1)*dtheta(2)*(1-cos(angle)))-(dtheta(3)*sin(angle)) (dtheta(1)*dtheta(3)*(1-cos(angle)))+(dtheta(2)*sin(angle));...
            (dtheta(2)*dtheta(1)*(1-cos(angle)))+(dtheta(3)*sin(angle)) cos(angle)+((dtheta(2)^2)*(1-cos(angle))) (dtheta(2)*dtheta(3)*(1-cos(angle)))-(dtheta(1)*sin(angle));...
            (dtheta(3)*dtheta(1)*(1-cos(angle)))-(dtheta(2)*sin(angle)) (dtheta(3)*dtheta(2)*(1-cos(angle)))+(dtheta(1)*sin(angle)) cos(angle)+((dtheta(3)^2)*(1-cos(angle)))];
        a(step+1,:) = (R*(a(step,:)-c(step,:))')+c(step,:)'+(c(step+1,:)-c(step,:))';%+a(step,:)'%+(c(step,:)-a(step,:));%cross(dtheta,c(step+1,:));%a(step,:) + (c(step+1,:)-c(step,:))+cross(dtheta/norm(dtheta),c(step,:));%+c(step,:);
        
        tempLenMot=sqrt((a(step+1,1)-x_m(step,1))^2+(a(step+1,2)-x_m(step,2))^2+(a(step+1,3)-x_m(step,3))^2);
        f_mot_head=-k_m*(tempLenMot-l_m)*(tempLenMot>l_m)*((x_m(step,:)-a(step+1,:))/norm(x_m(step,:)-a(step+1,:)));
        x_m_part = x_m(step,:) + (dt*(1/gamma_m)*(f_mot_head+f_tail_on_head));
        acceptMove = 0;
        ct = 1;
        while acceptMove == 0
            x_m(step+1,:) = x_m_part+((sqrt(2*D_m))*sqrt(dt)*randn(1,3)); %Wiener process for c terminal tail
            %x_m(step+1,:) = x_m(step,:) + (dt*(1/gamma_m)*(f_mot_head+f_tail_on_head))+((sqrt(2*D_m))*W_m);%(dt*(-(k_m/gamma_m)*(norm(x_m(step,:)-a)-l_m)*((x_m(step,:)-a)/norm(x_m(step,:)-a))+k_b*(x_t(boundTailInd(2),:,boundTailInd(1))-x_m(step,:)))) + ((sqrt(2*D_m))*W_m);
            if sqrt(x_m(step+1,2)^2+x_m(step+1,3)^2)>=12.5 && sqrt(((c(step+1,1)-x_m(step+1,1))^2)+((c(step+1,2)-x_m(step+1,2))^2)+((c(step+1,3)-x_m(step+1,3))^2))>=car_radi && x_m(step+1,2)>-3
                acceptMove = 1;
            else
                ct = ct+1;
                if ct>10000
                    break
                end
            end
        end
        if ct>10000
            disp('CTT-bound motor cannot move')
%             step = 1;
%             t = zeros(step_max,1);
%             continue
            t(step) = time_max-1;
            MTbound=1;
            break
        end
%         if movie == 1 && mod(step,capFreq)==0
%             clf; 
%             mt = surf(Xc,Yc,Zc);hold on;
%             rotate(mt,[0 1 0], 90)
%             for xt = 1:numCTT
%                 for yt=1:numP+1
%                     plot3([x_t(yt,1,xt) b(yt,1,xt)], [x_t(yt,2,xt) b(yt,2,xt)],[x_t(yt,3,xt) b(yt,3,xt)],'color',[0.8500 0.3250 0.0980],'linewidth',2);
%                 end
%             end
%             plot3(a(step+1,1),a(step+1,2),a(step+1,3),'b.','MarkerSize',20)
%             plot3([x_m(step+1,1) a(step+1,1)], [x_m(step+1,2) a(step+1,2)],[x_m(step+1,3) a(step+1,3)],'linewidth',2);
%             plot3([x_t(boundTailInd(2),1,boundTailInd(1)) x_m(step+1,1)],[x_t(boundTailInd(2),2,boundTailInd(1)) x_m(step+1,2)],[x_t(boundTailInd(2),3,boundTailInd(1)) x_m(step+1,3)],...
%                 'color',[0.3010 0.7450 0.9330],'linewidth',3)
%             %plot(Cargo_x, Cargo_y);
%             surf(X+c(step+1,1),Y+c(step+1,2),Z+c(step+1,3),'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none')
%             view(2);axis equal
%             xlim([-200 200]);ylim([-12.5 300]);
%             title(['t= ' num2str(t(step)) 's'])
%             drawnow;
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         end
    else
%         if x_m(step,2)<12.5+l_t+bindRadiusCTT && abs(x_m(step,3))<12.5+l_t+bindRadiusCTT
%             if round((x_m(step,1)-CTT_xStart)/8) > 2
%                 num_it_min = round((x_m(step,1)-CTT_xStart)/8)-2;
%             else
%                 num_it_min = 1;
%             end
%             if round((x_m(step,1)-CTT_xStart)/8)<=numCTT-2
%                 num_it_max = round((x_m(step,1)-CTT_xStart)/8)+2;
%             elseif round((x_m(step,1)-CTT_xStart)/8) == numCTT
%                 num_it_max = round((x_m(step,1)-CTT_xStart)/8)+1;
%             else
%                 num_it_max = round((x_m(step,1)-CTT_xStart)/8);
%             end
%         else
%             num_it_min = randi(numCTT);
%             num_it_max = num_it_min;
%         end
        
%         for wt = num_it_min:num_it_max
%             for xt = 1:numP+1
%                 Px_t_part = x_t(xt,:,wt) + (-(k_t/gamma_t)*dt*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt))));
%                 acceptMove = 0;
%                 ct = 1;
%                 while acceptMove == 0
%                     Px_t = Px_t_part + ((sqrt(2*D_t))*sqrt(dt)*randn(1,3)); %Wiener process for c terminal tail
%                     %Px_t = x_t(xt,:,wt) + (-(k_t/gamma_t)*dt*(norm(x_t(xt,:,wt)-b(xt,:,wt))-l_t)*((x_t(xt,:,wt)-b(xt,:,wt))/norm(x_t(xt,:,wt)-b(xt,:,wt)))) + ((sqrt(2*D_t))*W_t);
%                     if sqrt(Px_t(1,2)^2+Px_t(1,3)^2) >= 12.5 && Px_t(1,2)>-3
%                         acceptMove = 1;
%                         x_t(xt,:,wt)=Px_t;
%                     else
%                         ct = ct + 1;
%                         if ct > 1000
%                             break
%                         end
%                     end
%                 end
%                 if ct > 1000
%                     disp('unbound CTT cannot move')
%                     step = 1;
%                     t = zeros(step_max,1);
%                     continue
%                 end
%             end
%         end
        c_part = (k_laser.*(c(1,:)-c(step,:))*dt)+(f_mot_head*dt);
        accept = 0;
        ct = 1;
        while accept == 0
            c(step+1,:)=c(step,:)+(sqrt(dt*2*D_c)*randn(1,3))+ c_part;%(k_laser.*(c(1,:)-c(step,:))*dt)+(f_mot_head*dt);
            if c(step+1,2)>12.5+car_radi %&& abs(c(step+1,1))<abs(CTT_xStart)-car_radi-l_m-100
                accept = 1;
            else
                ct = ct + 1;
                if ct > 1000
                    break
                end
            end
        end
        if ct >1000
            disp('cargo in unbound case cannot move')
            step = 1;
            t = zeros(step_max,1);
            continue
        end
        if MTbound == 0
            torque_m = 0;
        else
            torque_m = cross((a(step,1)-c(step+1,:)),f_mot_head);
        end
        theta(step+1,:) = theta(step,:)+((1/(8*pi*eta*(car_radi^3)))*torque_m*dt)+(sqrt(2*(4.114/(8*pi*eta*(car_radi^3)))*dt)*randn(1,3));
        dtheta=theta(step+1,:)-theta(step,:);
        angle = norm(dtheta);
        tempAngle(step) = angle;
        dtheta = dtheta/angle;
        R = [cos(angle)+((dtheta(1)^2)*(1-cos(angle))) (dtheta(1)*dtheta(2)*(1-cos(angle)))-(dtheta(3)*sin(angle)) (dtheta(1)*dtheta(3)*(1-cos(angle)))+(dtheta(2)*sin(angle));...
            (dtheta(2)*dtheta(1)*(1-cos(angle)))+(dtheta(3)*sin(angle)) cos(angle)+((dtheta(2)^2)*(1-cos(angle))) (dtheta(2)*dtheta(3)*(1-cos(angle)))-(dtheta(1)*sin(angle));...
            (dtheta(3)*dtheta(1)*(1-cos(angle)))-(dtheta(2)*sin(angle)) (dtheta(3)*dtheta(2)*(1-cos(angle)))+(dtheta(1)*sin(angle)) cos(angle)+((dtheta(3)^2)*(1-cos(angle)))];
        a(step+1,:) = (R*(a(step,:)-c(step,:))')+c(step,:)'+(c(step+1,:)-c(step,:))';%+(c(step,:)-a(step,:));%cross(dtheta,c(step+1,:));%a(step,:) + (c(step+1,:)-c(step,:))+cross(dtheta/norm(dtheta),c(step,:));%+c(step,:);
        tempLenMot=sqrt((a(step+1,1)-x_m(step,1))^2+(a(step+1,2)-x_m(step,2))^2+(a(step+1,3)-x_m(step,3))^2);
        f_mot_head=-k_m*(tempLenMot-l_m)*(tempLenMot>l_m)*((x_m(step,:)-a(step+1,:))/norm(x_m(step,:)-a(step+1,:)));
        if MTbound == 1
            f_weak = k_weak*(weak_anchor-x_m(step,:));
        else
            f_weak = 0;
        end
        x_m_part = x_m(step,:) + (dt*(1/gamma_m)*(f_mot_head+f_weak));
        acceptMove = 0;
        ct = 1;
        while acceptMove == 0
            x_m(step+1,:) = x_m_part + ((sqrt(2*D_m))*sqrt(dt)*randn(1,3)); %Wiener process for motor
            %x_m(step+1,:) = x_m(step,:) + (dt*(1/gamma_m)*f_mot_head)+((sqrt(2*D_m))*W_m);
            if sqrt(x_m(step+1,2)^2+x_m(step+1,3)^2) >= 12.5 && sqrt(((c(step+1,1)-x_m(step+1,1))^2)+((c(step+1,2)-x_m(step+1,2))^2)+((c(step+1,3)-x_m(step+1,3))^2))>=car_radi && x_m(step+1,2)>-3
                acceptMove=1;
            else
                ct = ct+1;
                if ct>10000
                    break
                end
            end
        end
        if ct>10000
            disp('free motor cannot move')
%             step = 1;
%             t = zeros(step_max,1);
%             continue
            t(step) = time_max-1;
            MTbound=1;
            break
        end
        if movie == 1 && (mod(step,capFreq)==0 || step==1)
            clf; 
            mt = surf(Xc,Yc,Zc);hold on
            rotate(mt,[0 1 0], 90)
%             for xt = 1:numCTT
%                 for yt=1:numP+1
%                     plot3([x_t(yt,1,xt) b(yt,1,xt)], [x_t(yt,2,xt) b(yt,2,xt)],[x_t(yt,3,xt) b(yt,3,xt)],'color',[0.8500 0.3250 0.0980],'linewidth',2);
%                 end
%             end
            plot3(a(step,1),a(step,2),a(step,3),'b.','MarkerSize',20)
            plot3([x_m(step,1) a(step,1)], [x_m(step,2) a(step,2)],[x_m(step,3) a(step,3)],'linewidth',2);
            surf(X+c(step,1),Y+c(step,2),Z+c(step,3),'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none');
            view(2);axis equal
            xlim([-200 200]);ylim([-12.5 300]);
            title(['t = ' num2str(t(step)) 's'])
            %xlim([-a(2) a(2)]);ylim([-12.5 2*a(2)-12.5]);zlim([-a(2) a(2)])
            drawnow;
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
    
    closeTails = [];
%     for xt = num_it_min:num_it_max
%         for yt = 1:numP+1
%             distMotorCTT(xt,yt) = sqrt((x_m(step+1,1)-x_t(yt,1,xt))^2 + (x_m(step+1,2)-x_t(yt,2,xt))^2 + (x_m(step+1,3)-x_t(yt,3,xt))^2);
%             if distMotorCTT(xt,yt)<bindRadiusCTT
%                 closeTails=[closeTails; xt yt];
%             end
%         end
%     end
    LenMot=sqrt((a(step+1,1)-x_m(step+1,1))^2+(a(step+1,2)-x_m(step+1,2))^2+(a(step+1,3)-x_m(step+1,3))^2);
%     disp(f_mot_head)
%     disp(tailbound)
%     disp(LenMot)
%     disp(x_m)
    step = step + 1;
    t(step) = t(step-1) + dt;
    
    if t(step)>=time_max || step == step_max 
        %disp('Did not bind');
        step = 1;
        t = zeros(step_max,1);
        continue
    end
end
if movie==1
    close(v);
end
if t(step)<time_max && MTbound==1
    timetobind = t(step);
    tFinal = t(1:step);
    %disp(timeTailbound);
    %timeTailbound = timeTailbound/timetobind;
    %disp(t(step));
    %j = j+1
else
    disp('Did not bind');
    ct
end
% figure(2);clf;hold on;
% % % plot(t,distMotorCTT)
%histogram(MarchovState)
% MT=plot(linspace(0,120,1000),12.5*ones(1000,1),'k--');
% C=plot(linspace(0,120,1000),a(2)*ones(1000,1),'r--');
% ax = gca;
% ax.LineWidth = 1.5;
% set(ax, 'fontsize', 15', 'color', 'w')
% legend([MT,C],'Top of Microtubule','Bottom-most of Cargo');
% legend('location', 'best')
% h3.Color = 'w';
% title({['Distance between Motor and MT: ' num2str(motorTailSpace+l_t) ' nm,'];...
%     ['Motor-MT On-rate = ' num2str(k_on_MT) ' s^{-1}, Motor-CTT On-rate = ' num2str(k_on_CTT) ' s^{-1}']})
% ylabel('Time to Bind (sec)', 'color', 'k')
% xlabel('k^{ADP}_{off} (s^{-1})', 'color', 'k')
% box on;grid on

%end
 end

