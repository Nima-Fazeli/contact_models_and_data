
clc
clear
close all

addpath('../bag_files')

withVis  = 0;
withSave = 0;
second_missing = [];

start_num = 508;
numData   = 2000;

bounce_array(1).states = zeros(1,12);
bounce_array(1).n = zeros(1,3);
bounce_array(1).d = zeros(1,3);
bounce_array(1).flag_sb = 0;
bounce_array(1).flag_en = 0;
bounce_array(1).flag_cl = 0;
bounce_array(1).flag_mc = 0;
bounce_array(1).flag_de = 0;
bounce_array(1).flag_ld = 0;
bounce_array(1).flags = [0,0,0,0,0,0];
clean_data = 0;

for count = start_num:numData
    close all
    flag_sb = 0; % second bounce
    flag_en = 0; % energy
    flag_cl = 0; % second bounce
    flag_mc = 0; % second bounce
    flag_de = 0; % dropped early
    flag_ld = 0; % dropped frames, loss of data
    
    bag_num = count;
    bag_id = sprintf('drop_num_%d.bag',bag_num);
    bag = rosbag(bag_id);
    
    poseMsg = select(bag,'Topic','/object_arena');
    
    allMsg = readMessages(poseMsg);
    
    N = size(allMsg,1);
    
    p = zeros(N,3);
    q = zeros(N,4);
    r = zeros(N,3);
    
    for i=1:N
        p(i,:)=[allMsg{i}.Transform.Translation.X,...
            allMsg{i}.Transform.Translation.Y,...
            allMsg{i}.Transform.Translation.Z];
        
        q(i,:)=[allMsg{i}.Transform.Rotation.X,...
            allMsg{i}.Transform.Rotation.Y,...
            allMsg{i}.Transform.Rotation.Z,...
            allMsg{i}.Transform.Rotation.W];
        
        r(i,:) = quat2eul(q(i,:));
    end
    
    if withVis
        %%% plot raw data
        figure
        subplot(211)
        hold on
        for i=1:3
            plot(p(1:400,i))
        end
        legend('x','y','z','Location','NorthWest')
        grid
        box
        
        subplot(212)
        hold on
        for i=1:3
            plot(r(1:400,i))
        end
        legend('r_x','r_y','r_z','Location','NorthWest')
        grid
        box
    end
    
    %%% unwrap the angles
    threshold = 179;
    r_c=[unwrap(r(:,1),threshold*pi/180),...
        unwrap(r(:,2),threshold*pi/180),...
        unwrap(r(:,3),threshold*pi/180)];
    
    if withVis
        %%% plot raw data with unwrapped angles
        figure
        subplot(211)
        hold on
        for i=1:3
            plot(p(1:400,i))
        end
        legend('x','y','z','Location','NorthWest')
        grid
        box
        
        subplot(212)
        hold on
        for i=1:3
            plot(r_c(1:400,i))
        end
        legend('r_x','r_y','r_z','Location','NorthWest')
        grid
        box
    end
    
    %%% compute and plot velocities
    fs = 250;
    
    vl = diff(p)*fs;
    vr = diff(r_c)*fs;
    
    if withVis
        figure
        subplot(211)
        hold on
        for i=1:3
            plot(vl(1:400,i))
        end
        legend('x','y','z','Location','NorthWest')
        grid
        box
        
        subplot(212)
        hold on
        plot(vr(1:400,1))
        legend('r_x','Location','NorthWest')
        grid
        box
    end
    
    %%% compute total energy
    a = 70e-3/2;
    b = 50e-3/2;
    Im = (a^2+b^2)/5;
    g = 9.81;
    M = diag([1,1,Im]);
    
    dotq = [vl(:,2),vl(:,3),vr(:,1)];
    
    T = zeros(N-1,1);
    U = zeros(N-1,1);
    for i=1:N-1
        T(i) = 0.5*dotq(i,:)*M*dotq(i,:)';
        U(i) = g*p(i,3);
    end
    
    if withVis
        figure
        subplot(311)
        plot(T(1:400))
        box
        grid
        ylabel 'Kinetic'
        subplot(312)
        plot(U(1:400))
        box
        grid
        ylabel 'Potential'
        subplot(313)
        plot(U(1:400)+T(1:400))
        box
        grid
        ylabel 'Total'
    end
    
    E = U+T;
    
    %%%% investigate acceleration
    
    acc_threshhold = 200;
    Rg = sqrt(Im);
    acc_vec = [diff(p(:,2),2),diff(p(:,3),2),diff(r_c(:,1),2)*Rg]*fs*fs;
    acc_norm = sqrt(sum(abs(acc_vec).^2,2));
    if withVis
        figure
        subplot(311)
        plot(p(:,3))
        subplot(312)
        plot(diff(p(:,3))*fs)
        subplot(313)
        hold on
        plot(diff(p(:,3),2)*fs*fs)
        plot(acc_norm)
    end
    
    bounce_ind = find(acc_norm>200);
    if isempty(bounce_ind)
        bounce_ind = find(acc_norm>150);
        if isempty(bounce_ind)
            bounce_array(count).states = drop_states;
            bounce_array(count).n = n_for_storage';
            bounce_array(count).d = d_for_storage';
            bounce_array(count).flag_sb = 1;
            bounce_array(count).flag_en = 1;
            bounce_array(count).flag_cl = 1;
            bounce_array(count).flag_mc = 1;
            bounce_array(count).flag_de = 1;
            bounce_array(count).flag_ld = 1;
            bounce_array(count).flags = [1,1,1,1,1,1];
            continue
        end
    end
    bounce_ind = bounce_ind+1;
    
    if (bounce_ind(1)-18)>1
        pre_imp_mean = mean(dotq(bounce_ind(1)-5:bounce_ind(1)-1,2));
        pst_imp_mean = mean(dotq(bounce_ind(1)+2:bounce_ind(1)+5,2));
        if pre_imp_mean*pst_imp_mean < 0
            ind_beg = bounce_ind(1) - 20;
            if ind_beg < 1
                ind_beg = 1;
            end
            ind_end = bounce_ind(1) - 5;
        else
            flag_ld = 1;
            for iii=2:length(bounce_ind)
                if bounce_ind(iii)-bounce_ind(1)>7
                    ind_beg = bounce_ind(iii) - 20;
                    if ind_beg < 1
                        ind_beg = 1;
                    end
                    ind_end = bounce_ind(iii) - 5;
                    bounce_ind(1)=bounce_ind(iii);
                    break
                end
            end
        end
    else
        flag_de = 1;
        for find_first = 1:length(bounce_ind)
            if bounce_ind(find_first)>10
                ind_beg = bounce_ind(find_first) - 20;
                if ind_beg < 1
                    ind_beg = 1;
                end
                ind_end = bounce_ind(find_first) - 5;
                bounce_ind(1)=bounce_ind(find_first);
                break
            else
                continue
            end
        end
    end
    
    
    t = [0:1/250:(ind_end-ind_beg)/250]';
    yy = p(ind_beg:ind_end,3)+0.5*9.81*t.^2;
    yx = p(ind_beg:ind_end,2);
    yt = r_c(ind_beg:ind_end,1);
    x = [t,ones(length(t),1)];
    
    icy = lsqlin(x,yy);
    icx = lsqlin(x,yx);
    ict = lsqlin(x,yt);
    impy = [-0.5*9.81*(20/250)^2+icy(1)*(20/250)+icy(2);-9.81*(19/250)+icy(1)];
    impx = [icx(1)*(20/250)+icx(2);icx(1)];
    impt = [ict(1)*(20/250)+ict(2);ict(1)];
    
    if withVis
        subplot(311)
        hold on
        plot(ind_beg,icy(2),'or')
        plot(bounce_ind(1),impy(1),'ob')
        subplot(312)
        hold on
        plot(ind_beg,icy(1),'or')
        plot(bounce_ind(1)-1,impy(2),'ob')
    end
    
    %%% contact configuration
    ind_imp = bounce_ind(1)-1;
    pos = p(ind_imp,2:3);
    ori = r_c(ind_imp,1);
    ellipse = [70e-3,50e-3]./2;
    res = 1000;
    ic = draw_ellipse_func(  pos, ori, ellipse, res );
    [ Y_min, beta_c, X_min, contLoc] = min_point_ellipse( pos, ori, ellipse, res );
    
    conf = [pos';ori];
    n = calcN([X_min Y_min], conf);
    d = calcD([X_min Y_min], conf);
    
    n_for_storage = n;
    d_for_storage = d;
    
    J = [d';...
        n'];
    v_c = J*[vl(ind_imp-1,2:3)';vr(ind_imp-1,1)];
    m_c = inv(J*(M\J'));
    
    if withVis
        figure
        hold on
        plot(ic(1,:),ic(2,:))
        plot(X_min,Y_min,'or')
        axis equal
        grid
        box
    end
    
    %%% contact location mass matrix
    ind_pst = bounce_ind(1)+2;
    pos = p(ind_pst,2:3);
    ori = r_c(ind_pst,1);
    ellipse = [70e-3,50e-3]./2;
    res = 1000;
    ic = draw_ellipse_func(  pos, ori, ellipse, res );
    [ Y_min, beta_c, X_min, contLoc] = min_point_ellipse( pos, ori, ellipse, res );
    
    conf = [pos';ori];
    n = calcN([X_min Y_min], conf);
    d = calcD([X_min Y_min], conf);
    
    J = [d';...
        n'];
    
    v_p = J*[vl(ind_pst,2:3)';vr(ind_pst,1)];
    
    energy_check = v_c'*m_c*v_c - v_p'*m_c*v_p;
    
    rad = v_c'*m_c*v_c;
    MVI = m_c*v_c;
    
    [xe,ye,minmax] = drawE2(m_c, MVI, rad);
    
    xaxp=linspace(0,1.5,100);
    xaxn=linspace(-1.5,0,100);
    
    P = m_c*(v_p-v_c);
    
    mu = 0.35;
    
    if withVis
        figure
        hold on
        plot(xe,ye,'b')
        plot(xaxp,1/mu*xaxp,'r')
        plot(xaxn,-1/mu*xaxn,'r')
        plot(P(1),P(2),'ok')
        plot([xe(minmax(1)),xe(minmax(2))],[ye(minmax(1)),ye(minmax(2))],'k')
        box
        grid
        axis equal
    end
    
    %%% impact and the second parabola
    fprintf('\n Looking at trial number %d: \n',count)
    
    fb = bounce_ind(1);
    sb =[];
    
    if length(bounce_ind)<2
        fprintf('WARNING: I have detected 1 impact!!! \n')
        second_missing = [second_missing,count];
        sb = fb+5;
        flag_sb = 1;
    else
        fprintf('I have detected %d impacts. \n',length(bounce_ind))
        for ii=1:length(bounce_ind)
            fprintf('Impact %d is at %d\n',ii,bounce_ind(ii))
        end
    end
    
    for ii=2:length(bounce_ind)
        if bounce_ind(ii)-fb <7
            continue
        else
            sb = bounce_ind(ii);
            break
        end
    end
    
    if isempty(sb)
        fprintf('False positive second bounce, skipping \n');
        sb = fb+5;
        flag_sb = 1;
    end
    
    ind_para_a = fb+2;
    ind_para_b = sb-2;
    
    t = [0:1/250:(ind_para_b-ind_para_a)/250]';
    y = p(ind_para_a:ind_para_b,3)+0.5*9.81*t.^2;
    x = [t,ones(length(t),1)];
    
    ic = lsqlin(x,y);
    pimp = [-0.5*9.81*(-2/250)^2+ic(1)*(-2/250)+ic(2);-9.81*(-1/250)+ic(1)];
    
    if withVis
        figure
        subplot(312)
        hold on
        x_axis = ind_para_a:ind_para_b;
        t_axis = 0:1/250:(ind_para_b-ind_para_a)/250;
        plot(x_axis,-9.81*t_axis+ic(1),'--r')
        plot(ind_para_a-1,pimp(2),'or')
    end
    
    t = [0:1/250:(ind_para_b-ind_para_a)/250]';
    yy = p(ind_para_a:ind_para_b,3)+0.5*9.81*t.^2;
    yx = p(ind_para_a:ind_para_b,2);
    yt = r_c(ind_para_a:ind_para_b,1);
    x = [t,ones(length(t),1)];
    
    icy = lsqlin(x,yy);
    icx = lsqlin(x,yx);
    ict = lsqlin(x,yt);
    pimpy = [-0.5*9.81*(-2/250)^2+icy(1)*(-2/250)+ic(2);-9.81*(-1/250)+icy(1)];
    pimpx = [icx(1)*(-2/250)+icx(2);icx(1)];
    pimpt = [ict(1)*(-2/250)+ict(2);ict(1)];
    
    
    %%% data formating
    drop_states = [impx(1),impy(1),impt(1),impx(2),impy(2),impt(2),...
        pimpx(1),pimpy(1),pimpt(1),pimpx(2),pimpy(2),pimpt(2)];
    
    
    %%%% tests
    %%% flag 1, energy test
    if energy_check < 0
        fprintf('Energy test failed at %d\n',count);
        flag_en = 1;
    end
    
    %%% flag 2, coloumb test
    c1 =P(2)-1/mu*P(1);
    c2 =P(2)+1/mu*P(1);
    
    if c1 < 0 || c2 <0
        fprintf('Friction cone test failed at %d\n',count);
        flag_cl = 1;
    end
    
    %%% flag 3, maximum compression test
    m=(ye(minmax(2))-ye(minmax(1)))/(xe(minmax(2))-xe(minmax(1)));
    max_comp_cons = P(2)-ye(minmax(1))-m*(P(1)-xe(minmax(1)));
    
    if max_comp_cons < 0
        fprintf('Max compression constraint violated at %d\n',count);
        flag_mc = 1;
    end
    
    clean_data = clean_data + 1;
    bounce_array(count).states = drop_states;
    bounce_array(count).n = n_for_storage';
    bounce_array(count).d = d_for_storage';
    bounce_array(count).flag_sb = flag_sb;
    bounce_array(count).flag_en = flag_en;
    bounce_array(count).flag_cl = flag_cl;
    bounce_array(count).flag_mc = flag_mc;
    bounce_array(count).flag_de = flag_de;
    bounce_array(count).flag_ld = flag_ld;
    bounce_array(count).flags = [flag_sb,flag_en,flag_cl,flag_mc,flag_de,flag_ld];
    
    if withSave
        cd ../data
        save('ellipse_uniform','bounce_array')
        cd ../data_process
    end
    
    if rem(count,10)==0
        fprintf('\nPercent done: %2.2f',count/numData*100)
    end
end

