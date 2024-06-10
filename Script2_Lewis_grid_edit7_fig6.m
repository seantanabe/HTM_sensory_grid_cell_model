

clear all

N_sense_type            = 40;
N_sense_type_per_object = 10;
N_sense_col             = 150;
N_cell_per_col          = 16;
N_sense_per_object      = 10;
N_object_sq             = 4;
N_object                = 50; %100;

sense_sdr = {};
for i_sen = 1:N_sense_type
    sen_tmp = zeros(N_sense_col,1);
    sense_ind = randperm(length(sen_tmp(:)),10);
    sen_tmp(sense_ind) = 1;
    sense_sdr{i_sen} = sen_tmp;
end

objects = {};
for i_obj = 1:N_object
    obj_tmp = zeros(N_object_sq,N_object_sq);
    sense_ind = randperm(length(obj_tmp(:)),N_sense_per_object);
    sense_type_i = randperm(N_sense_type, N_sense_type_per_object);
    sense_type_ind = randi(N_sense_type_per_object,N_sense_per_object,1);
    obj_tmp(sense_ind) = sense_type_i(sense_type_ind);
    objects{i_obj} = obj_tmp;
end

N_grid_module = 10;
N_grid_sq = 40; %6;

object_width = 10;
scale = object_width/2;
object_sense_pos = 0:(object_width/N_object_sq):object_width;
sense_center_pos = (object_sense_pos(2)/2):object_sense_pos(2):object_width;
theta_module = {}; M_module = {};
for i_m = 1:N_grid_module
    theta_module{i_m} = pi*(1/3)*rand;
    M_module{i_m} = inv([scale*cos(theta_module{i_m}) scale*cos(theta_module{i_m} + deg2rad(60)); ...
        scale*sin(theta_module{i_m}) scale*sin(theta_module{i_m} + deg2rad(60))]);
end

sigma = 0.18172/(N_grid_sq*(1/6));
deltaTheta = (1/3)/(N_grid_sq*(1/6));
a_active = exp(-(((deltaTheta/2)*(2/sqrt(3)))^2)/(2*(sigma^2)));

thrsh_dend_loc   = 8;
thrsh_dend_sense = ceil((N_grid_module)*0.8);

grid_len = 1/N_grid_sq;
phase_c = {}; % figure
for i_c = 1:N_grid_sq
    for ii_c = 1:N_grid_sq
        grid_len_i = grid_len/2 + grid_len*(i_c-1);
        grid_len_ii = grid_len/2 + grid_len*(ii_c-1);
        phase_c{i_c,ii_c} = [grid_len_i grid_len_ii];

        %         scatter(phase_c{i_c,ii_c}(1),phase_c{i_c,ii_c}(2),'.k'); hold on
    end
end

rand_pos = 0;
if rand_pos == 1
    N_step = 100; %1000; %5000;
    % step_max = object_width/N_grid_sq;
    step_max = object_width/N_grid_sq + 1;
    pos_mouse = zeros(2,N_step); pos_ini = [object_width/2 object_width/2];
    for i_st = 1:N_step
        step_i = rand([2,1])*2*step_max - step_max;
        if i_st == 1
            pos_mouse(:,i_st) = pos_ini;
        else
            pos_mouse(:,i_st) = pos_mouse(:,i_st-1) + step_i;
        end
        if pos_mouse(1,i_st) > object_width | pos_mouse(1,i_st) < 0
            pos_mouse(1,i_st) = pos_mouse(1,i_st-1) - step_i(1);
        end
        if pos_mouse(2,i_st) > object_width | pos_mouse(2,i_st) < 0
            pos_mouse(2,i_st) = pos_mouse(2,i_st-1) - step_i(2);
        end
    end
end

sense_pos = 1;
if sense_pos == 1
    pos_incr = object_sense_pos(2);
    pos_mouse = zeros(2,N_object_sq^2);
    count = 0;
    step_i = repelem(1:N_object_sq,4);
    step_ii = [1:N_object_sq N_object_sq:-1:1 1:N_object_sq N_object_sq:-1:1];
    for i_step = 1:(N_object_sq^2)
        count = count + 1;
        pos_mouse(:,count) = [(object_sense_pos(2)/2+(step_i(i_step)-1)*pos_incr) (object_sense_pos(2)/2+(step_ii(i_step)-1)*pos_incr)];
    end
    pos_mouse = [pos_mouse(:,2) pos_mouse];
end

% figure; plot(pos_mouse(1,1:500),pos_mouse(2,1:500)); xlim([0 object_width]); ylim([0 object_width])
% pos_mouse = rand([2,10000])*object_width; % does not work, maybe mod1 lim
d_mouse   = diff(pos_mouse,[],2);

D_sense_loc = zeros(N_sense_col*N_cell_per_col, N_grid_module*(N_grid_sq^2))';
A_sense_rs_all = zeros(N_sense_col*N_cell_per_col, size(d_mouse,2)*length(objects));
A_loc_rs_all = zeros(N_grid_module*(N_grid_sq^2),size(d_mouse,2)*length(objects));
A_loc_all = zeros(N_grid_sq,N_grid_sq,N_grid_module,size(d_mouse,2)*length(objects));
for i_o = 1:length(objects)
    object_i      = objects{i_o};
    phase_loc_ini = rand(2,N_grid_module);
    for i_d = 1:size(d_mouse,2)
        % i_d = 1
        if rem(i_d,10) == 0
            disp(['object ' num2str(i_o) ' move ' num2str(i_d)])
        end
        if i_d == 1
            phase_loc = phase_loc_ini;
        end
        d = d_mouse([1 2],i_d)';

        a_cb  = zeros(N_grid_sq,N_grid_sq,N_grid_module);
        A_loc = zeros(N_grid_sq,N_grid_sq,N_grid_module);
        for i_m = 1:N_grid_module
            phase_loc(:,i_m) = mod(phase_loc(:,i_m) + M_module{i_m}*d',1);
            for i_c = 1:N_grid_sq
                for ii_c = 1:N_grid_sq
                    d_tmp = pdist([phase_c{i_c,ii_c}; phase_loc(:,i_m)'],'euclidean');
                    a_cb(i_c,ii_c,i_m)  = exp(-(d_tmp^2)/(2*(sigma^2)));
                    A_loc(i_c,ii_c,i_m) = a_cb(i_c,ii_c,i_m) >= a_active;
                end
            end
        end
        A_loc_rs = reshape(A_loc,[1 length(A_loc(:))]);  %%%%%%%%%%
        potent_sense_rs = A_loc_rs*D_sense_loc >= thrsh_dend_sense; %%%%%%%%%%%%%%%%
        potent_sense = reshape(potent_sense_rs, [N_sense_col N_cell_per_col]); %%%%%%%%%%%%%%
        %
        %         %%% check reshape order
        %         A_loc_mod_id = zeros(size(A_loc));
        %         A_loc_grid_id = zeros(size(A_loc));
        %         for i_m = 1:N_grid_module;
        %             A_loc_mod_id(:,:,i_m) = i_m;
        %             A_loc_grid_id(:,:,i_m) = reshape(1:(N_grid_sq^2),N_grid_sq,N_grid_sq);
        %         end
        %         A_loc_mod_id_rs = reshape(A_loc_mod_id,[1 length(A_loc(:))]);
        %         A_loc_grid_id_rs = reshape(A_loc_grid_id,[1 length(A_loc(:))]);


        pos_i   = pos_mouse(:,i_d);

        obj_sen_pos_x = find((object_sense_pos < pos_i(1))); obj_sen_pos_x = obj_sen_pos_x(end);
        obj_sen_pos_y = find((object_sense_pos < pos_i(2))); obj_sen_pos_y = obj_sen_pos_y(end);
        sense_i = object_i(obj_sen_pos_y, obj_sen_pos_x);

        if sense_i == 0
            A_sense_col = zeros(N_sense_col,1);
        else
            A_sense_col = sense_sdr{sense_i};
        end
        A_sense        = zeros(N_sense_col,N_cell_per_col);
        A_sense_tmp    = zeros(N_sense_col,N_cell_per_col);
        for i_s = find(A_sense_col)'
            % i_s = find(A_sense_col); i_s = i_s(1);
            if ~any(potent_sense(i_s,:))
                A_sense(i_s,:) = 1;
                A_sense_tmp(i_s,randperm(N_cell_per_col,1)) = 1;
            else
                A_sense(i_s,:)     = potent_sense(i_s,:);
                A_sense_tmp(i_s,:) = potent_sense(i_s,:);
            end
        end
        A_sense_rs     = reshape(A_sense',[1 length(A_sense(:))]); %%%%%%%%%%%%%%%
        A_sense_tmp_rs = reshape(A_sense_tmp',[1 length(A_sense_tmp(:))]); %%%%%%%%%%%

        %%% check reshape order
        %     A_sense_col_id = zeros(size(A_sense));
        %     for i_s = 1:length(A_sense_col)
        %         A_sense_col_id(i_s,:) = i_s;
        %     end
        %     A_sense_col_id_rs = reshape(A_sense_col_id',[1 length(A_sense(:))]);

        D_sense_loc = max(D_sense_loc, A_sense_tmp_rs.*A_loc_rs'); % find(A_loc_rs) % find(A_sense_tmp_rs) %%%%%%%%%%%%%%%

        A_loc_all(:, :,:, i_d+size(d_mouse,2)*(i_o-1)) = A_loc;
        A_loc_rs_all(:, i_d+size(d_mouse,2)*(i_o-1)) = A_loc_rs;
        A_sense_rs_all(:, i_d+size(d_mouse,2)*(i_o-1)) = A_sense_rs;

    end
end

figure;
blur_sm = 5;
c_scale = 0.3;
tmp_img = imgaussfilt(A_loc_rs_all,blur_sm );
c_mxx = max(tmp_img(:))*c_scale;
subplot(2,3,1); imagesc(tmp_img); clim([0 c_mxx]); colormap(flip(gray)); xlabel('step'); ylabel('grid cell (sorted module)')
title('Active grid cells')
tmp_img = imgaussfilt(A_sense_rs_all,blur_sm );
c_mxx = max(tmp_img(:))*c_scale;
subplot(2,3,2); imagesc(tmp_img); clim([0 c_mxx]); colormap(flip(gray)); xlabel('step'); ylabel('sensor cell (sorted column)')
title('Active sensor cells')
tmp_img = imgaussfilt(D_sense_loc,blur_sm );
c_mxx = max(tmp_img(:))*c_scale;
subplot(2,3,3); imagesc(tmp_img); clim([0 c_mxx]); colormap(flip(gray)); ylabel('grid cell (sorted module)'); xlabel('sensor cell (sorted column)')
title('Synapses of dendrite')
subplot(2,3,4)
plot(pos_mouse(1,:),pos_mouse(2,:),'k'); hold on
scatter(pos_mouse(1,:),pos_mouse(2,:),50,'.r'); hold off; xlim([0 object_width]); ylim([0 object_width])
title('Mouse position')
set(gca,'YTick',[]); set(gca,'XTick',[])
subplot(2,3,5); imagesc(object_i);
set(gcf,'color','w');
title('Example object')
set(gca,'YTick',[]); set(gca,'XTick',[])

for mod_ind = 4:6 %N_grid_module;
    figure
    for step_i = 1:15 %20
        subplot(4,5,step_i); imagesc(A_loc_all(:, :,mod_ind, step_i)); clim([0 1]); colormap(flip(gray))
        title(['sense ' num2str(step_i)])
    end
    sgtitle(['module ' num2str(mod_ind)])
end
set(gcf,'color','w');
% 
% for mod_ind = 4:6 %N_grid_module;
%     figure
%     for step_i = 16:30 %20
%         subplot(4,5,step_i-15); imagesc(A_loc_all(:, :,mod_ind, step_i)); clim([0 1]); colormap(flip(gray))
%         title(['sense ' num2str(step_i)])
%     end
%     sgtitle(['2 module ' num2str(mod_ind)])
% end

%% inference

sense_rand_pos = 0;
if sense_rand_pos == 1
    N_step = 10;
    step_max = object_sense_pos(2);
    pos_mouse = zeros(2,N_step); pos_ini = [step_max/2 step_max/2];
    for i_st = 1:N_step
        step_i = randi([-1,1],2,1)*step_max;
        if i_st == 1
            pos_mouse(:,i_st) = pos_ini;
        else
            pos_mouse(:,i_st) = pos_mouse(:,i_st-1) + step_i;
        end
        if pos_mouse(1,i_st) > object_width | pos_mouse(1,i_st) < 0
            pos_mouse(1,i_st) = pos_mouse(1,i_st-1) - step_i(1);
        end
        if pos_mouse(2,i_st) > object_width | pos_mouse(2,i_st) < 0
            pos_mouse(2,i_st) = pos_mouse(2,i_st-1) - step_i(2);
        end
    end
end

sense_rand_pos_no0 = 1;
if sense_rand_pos_no0 == 1
    N_repeat = 1; 
    pos_mouse_no0 = {};
    for i_o = 1:N_object
        object_i      = objects{i_o};
        [ind_y, ind_x] = find(object_i > 0);
        pos_order_i = [];
        for i_r = 1:N_repeat
            pos_order_i = [pos_order_i randperm(N_sense_per_object)];
        end
        pos_mouse_no0{i_o} = [sense_center_pos(ind_x(pos_order_i)); sense_center_pos(ind_y(pos_order_i))];
    end
end

A_sense_rs_all = zeros(N_sense_col*N_cell_per_col, size(d_mouse,2)*length(objects));
A_loc_rs_all = zeros(N_grid_module*(N_grid_sq^2),size(d_mouse,2)*length(objects));
A_loc_all = zeros(N_grid_sq,N_grid_sq,N_grid_module,size(d_mouse,2)*length(objects));
for i_o = 1 %1:length(objects)
    object_i      = objects{i_o};
    pos_mouse = pos_mouse_no0{i_o} ;
    d_mouse   = diff(pos_mouse,[],2);
    phase_loc_ini = rand(2,N_grid_module);
    for i_d = 1:size(d_mouse,2)
        % i_d = 1
        if rem(i_d,1) == 0
            disp(['object ' num2str(i_o) ' move ' num2str(i_d)])
        end
        if i_d == 1
            phase_loc = phase_loc_ini;
        end
        d = d_mouse([1 2],i_d)';

        if i_d == 1
%             a_cb  = zeros(N_grid_sq,N_grid_sq,N_grid_module);
            A_loc = zeros(N_grid_sq,N_grid_sq,N_grid_module);
%             for i_m = 1:N_grid_module
%                 phase_loc(:,i_m) = mod(phase_loc(:,i_m) + M_module{i_m}*d',1);
%                 for i_c = 1:N_grid_sq
%                     for ii_c = 1:N_grid_sq
%                         d_tmp = pdist([phase_c{i_c,ii_c}; phase_loc(:,i_m)'],'euclidean');
%                         a_cb(i_c,ii_c,i_m)  = exp(-(d_tmp^2)/(2*(sigma^2)));
%                         A_loc(i_c,ii_c,i_m) = a_cb(i_c,ii_c,i_m) >= a_active;
%                     end
%                 end
%             end
        else
            A_loc = zeros(N_grid_sq,N_grid_sq,N_grid_module);
            for i_m = 1:N_grid_module
                potent_loc_i = potent_loc(:,:,i_m);
                [x_ii, y_ii] = find(potent_loc_i);
                phase_loc_ind = [x_ii y_ii];

                a_cb_cat  = zeros(N_grid_sq,N_grid_sq,size(phase_loc_ind,1));
                for i_p = 1:size(phase_loc_ind,1)
                    a_cb_i  = zeros(N_grid_sq,N_grid_sq);
                    phase_c_p = phase_c{phase_loc_ind(i_p,1),phase_loc_ind(i_p,2)}';
                    phase_c_p = mod(phase_c_p + M_module{i_m}*d',1);
                    for i_c = 1:N_grid_sq
                        for ii_c = 1:N_grid_sq
                            d_tmp = pdist([phase_c{i_c,ii_c}; phase_c_p'],'euclidean');
                            a_cb_i(i_c,ii_c)  = exp(-(d_tmp^2)/(2*(sigma^2)));
                        end
                    end
                    a_cb_cat(:,:,i_p) = a_cb_i;
                end
                a_cb = 1 - prod(1 - a_cb_cat,3) ;
                A_loc(:,:,i_m) = a_cb  >= a_active;
                %                 figure; imagesc(a_cb_cat_tmp(:,:,1))
            end
        end
        A_loc_rs = reshape(A_loc,[1 length(A_loc(:))]);  %%%%%%%%%%
        potent_sense_rs = A_loc_rs*D_sense_loc >= thrsh_dend_sense; %%%%%%%%%%%%%%%%
        potent_sense = reshape(potent_sense_rs, [N_sense_col N_cell_per_col]); %%%%%%%%%%%%%%
        %
        %         %%% check reshape order
        %         A_loc_mod_id = zeros(size(A_loc));
        %         A_loc_grid_id = zeros(size(A_loc));
        %         for i_m = 1:N_grid_module;
        %             A_loc_mod_id(:,:,i_m) = i_m;
        %             A_loc_grid_id(:,:,i_m) = reshape(1:(N_grid_sq^2),N_grid_sq,N_grid_sq);
        %         end
        %         A_loc_mod_id_rs = reshape(A_loc_mod_id,[1 length(A_loc(:))]);
        %         A_loc_grid_id_rs = reshape(A_loc_grid_id,[1 length(A_loc(:))]);


        pos_i   = pos_mouse(:,i_d);

        obj_sen_pos_x = find((object_sense_pos < pos_i(1))); obj_sen_pos_x = obj_sen_pos_x(end);
        obj_sen_pos_y = find((object_sense_pos < pos_i(2))); obj_sen_pos_y = obj_sen_pos_y(end);
        sense_i = object_i(obj_sen_pos_y, obj_sen_pos_x);

        if sense_i == 0
            A_sense_col = zeros(N_sense_col,1);
        else
            A_sense_col = sense_sdr{sense_i};
        end
        A_sense        = zeros(N_sense_col,N_cell_per_col);
        A_sense_tmp    = zeros(N_sense_col,N_cell_per_col);
        for i_s = find(A_sense_col)'
            % i_s = find(A_sense_col); i_s = i_s(1);
            if ~any(potent_sense(i_s,:))
                A_sense(i_s,:) = 1;
                A_sense_tmp(i_s,randperm(N_cell_per_col,1)) = 1;
            else
                A_sense(i_s,:)     = potent_sense(i_s,:);
                A_sense_tmp(i_s,:) = potent_sense(i_s,:);
            end
        end
        A_sense_rs     = reshape(A_sense',[1 length(A_sense(:))]); %%%%%%%%%%%%%%%
        A_sense_tmp_rs = reshape(A_sense_tmp',[1 length(A_sense_tmp(:))]); %%%%%%%%%%%

        %%% check reshape order
        %     A_sense_col_id = zeros(size(A_sense));
        %     for i_s = 1:length(A_sense_col)
        %         A_sense_col_id(i_s,:) = i_s;
        %     end
        %     A_sense_col_id_rs = reshape(A_sense_col_id',[1 length(A_sense(:))]);
        %
        % learning
        %         D_sense_loc = max(D_sense_loc, A_sense_tmp_rs.*A_loc_rs'); % find(A_loc_rs) % find(A_sense_tmp_rs) %%%%%%%%%%%%%%%

        %         potent_loc_rs = (D_sense_loc*A_sense_tmp_rs')' >= thrsh_dend_loc; %%%%%%%%%%%%%%%%
        potent_loc_rs = (D_sense_loc*A_sense_rs')' >= thrsh_dend_loc; %%%%%%%%%%%%%%%%
        potent_loc = reshape(potent_loc_rs, [N_grid_sq N_grid_sq N_grid_module]); %%%%%%%%%%%%%%
        %
        A_loc_all(:, :,:, i_d+size(d_mouse,2)*(i_o-1)) = A_loc;
        A_loc_rs_all(:, i_d+size(d_mouse,2)*(i_o-1)) = A_loc_rs;
        A_sense_rs_all(:, i_d+size(d_mouse,2)*(i_o-1)) = A_sense_rs;

    end
end

figure;
subplot(2,3,1); imagesc(A_loc_rs_all);   colormap(flip(gray)); xlabel('step'); ylabel('grid cell (sorted module)')
subplot(2,3,2); imagesc(A_sense_rs_all); colormap(flip(gray)); xlabel('step'); ylabel('sensor cell (sorted column)')
subplot(2,3,3); imagesc(D_sense_loc);    colormap(flip(gray)); ylabel('grid cell (sorted module)'); xlabel('sensor cell (sorted column)')
subplot(2,3,4)
plot(pos_mouse(1,:),pos_mouse(2,:),'k'); hold on
scatter(pos_mouse(1,:),pos_mouse(2,:),50,'.r'); hold off; xlim([0 object_width]); ylim([0 object_width])
subplot(2,3,5); imagesc(object_i);
set(gca,'YDir','normal')

for mod_ind = 4:6 %N_grid_module;
    figure
    for step_i = 1:10
        subplot(4,5,step_i); imagesc(A_loc_all(:, :,mod_ind, step_i)); clim([0 1]); colormap(flip(gray))
        title(['sense ' num2str(step_i)])
    end
    sgtitle(['module ' num2str(mod_ind)])
end

