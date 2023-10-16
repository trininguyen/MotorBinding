function [timetobind,MarchovState] = DiffwADPState(dist,k_ADP_off,k_ADP_on,k_ADP_fast_off,...
    k_MT_dir_on,k_MT_dir_off,k_CTT_on,k_CTT_off,k_CTT_MT,k_MT_CTT,length,D_m,k_weak,car_radi)

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
dtmax = 0.004; %difference of time between time steps (s)
t = zeros(step_max,1); %time

eta = 3e-5; %(pN/nm^2)*s
k_m = 0.320; %stiffness of kinesin (pN/nm)
k_laser = [0.045 0.03 0.045];%0.065; %stiffness of laser (pN/nm)
D_c = 4.114/(6*pi*eta*car_radi); %diffusion constant of cargo 
gamma_m=4.114/D_m;%motor drag coefficient

l_m = length+10; %length of kinesin motor (nm), +10 for binding antibody between motor and cargo (33.1429nm = 406aa, 45.7143nm=560aa, 60.8980=746aa)

bindRadiusMT = 12.5+2; %radius at which motor head can bind to MT (nm)

a = [0 12.5+dist 0]; %coordinate of motor anchor
x_m_init = [0 a(2)-l_m 0];
if x_m_init(2)<12.5
    x_m_init(2) = 12.5;
end

tempAngle = zeros(step_max,1);
theta = a;
c = [0 a(2)+car_radi 0]; %coordinate of cargo center
LenMot=sqrt((a(1)-x_m_init(1))^2+(a(2)-x_m_init(2))^2+(a(3)-x_m_init(3))^2);%length of motor
tempLenMot=LenMot;

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
    if step == 1
        tempTB = 0;
        tempTU = 0;
        trial = trial + 1;
        MotorADPState = 1; %if Motor is in ADP state, set to 1, if in ATP state, set to 2
        a = [0 12.5+dist 0]; %coordinate of motor anchor
        c = [0 a(2)+car_radi 0]; %coordinate of cargo center
        MTbound = 0; %if motor is bound to MT, set to 1
        f_mot_head = [0 0 0];
        x_m = zeros(step_max,3);
        x_m(1,:) = x_m_init;
        LenMot=sqrt((a(1)-x_m_init(1))^2+(a(2)-x_m_init(2))^2+(a(3)-x_m_init(3))^2);%length of motor
        tempLenMot=LenMot;
        MarchovState = [];
    end
    %Gillespie   
    if MTbound == 0
        if MotorADPState == 1
            ADP_off = k_ADP_off;
            ADP_on = 0;
            MT_on = k_MT_dir_on;
            MT_off = 0;
            MarchovState = [MarchovState 1];
        else
            ADP_off = 0;
            ADP_on = k_ADP_on;
            MT_on = k_MT_dir_on;
            MT_off = 0;
            MarchovState = [MarchovState 2];
        end
       if sqrt(x_m(step,2)^2+x_m(step,3)^2)>bindRadiusMT
           MT_on = 0;
       end
    else
        ADP_off = k_ADP_fast_off;
        ADP_on = 0;
        MT_on = k_ADP_fast_off;
        MT_off = k_MT_dir_off * exp(norm(f_weak)/4);
        MarchovState = [MarchovState 5];
        if MotorADPState==2
            ADP_off = 0;
        end
    end
    array_rates = [ADP_off ADP_on MT_on MT_off];
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
                MTbound = 1;
                if MotorADPState == 2 
                    if movie == 1
                        clf;
                        mt = surf(Xc,Yc,Zc);hold on;
                        rotate(mt,[0 1 0], 90)
                        plot3(a(step,1),a(step,2),a(step,3),'b.','MarkerSize',20)
                        plot3([x_m(step,1) a(step,1)], [x_m(step,2) a(step,2)],[x_m(step,3) a(step,3)],'linewidth',2);
                        surf(X+c(step,1),Y+c(step,2),Z+c(step,3),'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none')
                        view(2);axis equal
                        xlim([-200 200]);ylim([-12.5 300]);
                        title(['t = ' num2str(t(step)) 's'])
                        drawnow;
                        frame = getframe(gcf);
                        writeVideo(v,frame);
                    end
                    MarchovState = [MarchovState 7];
                    step = step + 1;
                    t(step) = t(step-1) + dt;
                    break;
                else
                   P = [x_m(step,2)-0,x_m(step,3)-0];
                    Q = (12.5/norm(P))*P;
                    weak_anchor = [x_m(step,1), Q(1) + x_m(step,2), Q(2) + x_m(step,3)];
                end
            case 4
                MTbound = 0;
                MotorADPState = 1;
        end
    else
       dt = dtmax;
    end
    
    c_part = (k_laser.*(c(1,:)-c(step,:))*dt)+(f_mot_head*dt);
    accept = 0;
    ct = 1;
    while accept == 0
        c(step+1,:)=c(step,:)+(sqrt(dt*2*D_c)*randn(1,3))+ c_part;
        if c(step+1,2)>12.5+car_radi 
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
    a(step+1,:) = (R*(a(step,:)-c(step,:))')+c(step,:)'+(c(step+1,:)-c(step,:))';
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
        t(step) = time_max-1;
        MTbound=1;
        break
    end
    if movie == 1 && (mod(step,capFreq)==0 || step==1)
        clf; 
        mt = surf(Xc,Yc,Zc);hold on
        rotate(mt,[0 1 0], 90)
        plot3(a(step,1),a(step,2),a(step,3),'b.','MarkerSize',20)
        plot3([x_m(step,1) a(step,1)], [x_m(step,2) a(step,2)],[x_m(step,3) a(step,3)],'linewidth',2);
        surf(X+c(step,1),Y+c(step,2),Z+c(step,3),'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none');
        view(2);axis equal
        xlim([-200 200]);ylim([-12.5 300]);
        title(['t = ' num2str(t(step)) 's'])
        drawnow;
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    LenMot=sqrt((a(step+1,1)-x_m(step+1,1))^2+(a(step+1,2)-x_m(step+1,2))^2+(a(step+1,3)-x_m(step+1,3))^2);
    step = step + 1;
    t(step) = t(step-1) + dt;
    
    if t(step)>=time_max || step == step_max 
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
else
    disp('Did not bind');
    ct
end
 end

