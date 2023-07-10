fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                      % 设置配置文件所在路径
% Verification_path = strcat(Config_Path,'sf10_BW500.txt');       % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf10_BW125.json'));    % 配置文件
Setting_File_Path = strcat(Config_Path, Setting_File.name);
Setting_file = fopen(Setting_File_Path,'r');
setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.bw
lora_set.sf = setting.captures.lora_sf;                         % 设置接收数据包的lora_set.sf
lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
% lora_set.sample_rate = 2e6;
lora_set.sample_rate = 500000;
lora_set.Pkg_length = setting.captures.lora_pkg_length;         % 设置接收数据包的长度
lora_set.dine = lora_set.sample_rate*bitshift(1,lora_set.sf)/lora_set.bw;    % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
lora_set.fft_x = 2^lora_set.sf;                               % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
lora_set.Preamble_length = 8;
leakage_width_array = [0.05,0.01,0.015,0.01];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;
lora_set.channel = 4;
lora_set.channel_choice = [0,0,0,0];
lora_set.channel_choice_num = 4;
lora_set.Preamble_channel = 0;

% 生成理想仿真信号 sf10,250KHz
fft_x = lora_set.fft_x;
dine = lora_set.dine;
lora_set.Preamble_channel = 0;
[~, upchirp_c0] = build_idealchirp_pre(lora_set);
lora_set.Preamble_channel = 1;
[~,upchirp_c1] = build_idealchirp_pre(lora_set);
lora_set.Preamble_channel = 2;
[~,upchirp_c2] = build_idealchirp_pre(lora_set);
lora_set.Preamble_channel = 3;
[~,upchirp_c3] = build_idealchirp_pre(lora_set);
% chirp = [upchirp_c0, upchirp_c1, upchirp_c2, upchirp_c3];
% ax1 = figure(1);
% stft(chirp, lora_set.sample_rate, 'Window',rectwin(128),'OverlapLength',127,'FFTLength',2048);
% colormap(ax1,CustomColormap);

% figure(1);
% cor = abs(xcorr(chirp, upchirp_c0, dine*4));
% cor = cor(dine*4+1:dine*8);
% plot(cor,'.r','LineWidth',6); hold on;
% cor = abs(xcorr(chirp, upchirp_c1, dine*4));
% cor = cor(dine*4+1:dine*8);
% plot(cor,'.b','LineWidth',6); hold on;
% cor = abs(xcorr(chirp, upchirp_c2, dine*4));
% cor = cor(dine*4+1:dine*8);
% plot(cor,'.g','LineWidth',6); hold on;
% cor = abs(xcorr(chirp, upchirp_c3, dine*4));
% cor = cor(dine*4+1:dine*8);
% plot(cor,'.m','LineWidth',6); hold on;

upchirp_array = zeros(4, dine);
dir_path = 'D:\CHchirp_samples\align_windows\rand_binchannel\channel4\preamble';
for i = 0:3
    file_path_tmp = strcat(dir_path, string(i),'\');
    file_name = dir(fullfile(file_path_tmp,'*.sigmf-data')); 
    file_path = strcat(file_path_tmp, file_name(1).name);
    fid_1=fopen(file_path,'rb');
    [A_1]=fread(fid_1,'float32')';
    A_1_length = size(A_1,2);
    G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;

    upchirp_array(i+1,:) = G0(1:4:dine*4-3);
end
chirp = [zeros(1,dine/2) ,upchirp_array(4,:), upchirp_array(3,:), upchirp_array(2,:), upchirp_array(1,:)];
% chirp = circshift(chirp, 30);
% ax1 = figure(1);
% stft(chirp, lora_set.sample_rate, 'Window',rectwin(128),'OverlapLength',127,'FFTLength',2048);
% colormap(ax1,CustomColormap);

figure(1);
cor = abs(xcorr(chirp, upchirp_c0, dine*4.5));
cor = cor(dine*4.5+1:dine*9);
plot(cor,'-','LineWidth',2,'color',[0 1 1]); hold on;
cor = abs(xcorr(chirp, upchirp_c1, dine*4.5));
cor = cor(dine*4.5+1:dine*9);
plot(cor,'--','LineWidth',2,'color',[0.07 0.62 1]); hold on;
cor = abs(xcorr(chirp, upchirp_c2, dine*4.5));
cor = cor(dine*4.5+1:dine*9);
plot(cor,'-.','LineWidth',2,'color',[0 0.45 0.74]); hold on;
cor = abs(xcorr(chirp, upchirp_c3, dine*4.5));
cor = cor(dine*4.5+1:dine*9);
plot(cor,':','LineWidth',2,'color',[0 0 1]); hold on;
ylabel('Correlation');
xlabel('PHY Samples');
xlim([0 dine*4.5]);
ylim([0 4500]);
set(gca,'FontSize',80);
set(get(gca,'YLabel'),'FontSize',100);
leg = legend('CH3','CH2','CH1','CH0');
leg.ItemTokenSize = [80,40];