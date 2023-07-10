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
lora_set.channel_choice = [3,2,1,0];
[d_upchirp_128, ~] = build_idealchirp_tmp(lora_set, 128);
lora_set.channel_choice = [3,1,2,0];
[d_upchirp_600, ~] = build_idealchirp_tmp(lora_set, 600);
lora_set.channel_choice = [0,1,2,3];
[d_upchirp_300, ~] = build_idealchirp_tmp(lora_set, 300);
upchirp_1 = [d_upchirp_128, d_upchirp_600];
upchirp_2 = [zeros(1,dine/2+800), d_upchirp_300, zeros(1,dine/2-800)];
dechirp = upchirp_1*0.5+upchirp_2;

% ax1 = figure(1);
% stft(dechirp(dine/2+800+1 : dine*3/2+800), lora_set.sample_rate, 'Window',rectwin(128),'OverlapLength',127,'FFTLength',2048);
% colormap(ax1,CustomColormap);
% colorbar('off');
% xlabel('Times(ms)');
% ylabel('Frequency(KHz)'); 
% ylim([-250 250]);
% yticks([-250 -125 0 125 250]);
% yticklabels({-250 -125 0 125 250});
% set(gca,'FontSize',80);
% set(get(gca,'YLabel'),'FontSize',100);
% set(gca,'xtick',[],'xticklabel',[]);


% ax2 = figure(1);
lora_set.channel_choice = [0,1,2,3];
[~, d_downchirp] = build_idealchirp_tmp(lora_set, 0);
dechirp = dechirp(dine/2+800+1 : dine*3/2+800).*d_downchirp;
% stft(dechirp, lora_set.sample_rate, 'Window',rectwin(256),'OverlapLength',0,'FFTLength',2048);
% colormap(ax2,CustomColormap);
% colorbar('off');
% xlabel('Times(ms)');
% ylabel('Frequency(KHz)'); 
% ylim([-250 250]);
% yticks([-250 -125 0 125 250]);
% yticklabels({-250 -125 0 125 250});
% set(gca,'FontSize',80);
% set(get(gca,'YLabel'),'FontSize',100);
% set(gca,'xtick',[],'xticklabel',[]);

% ax3 = figure(1);
% dechirp_fft = abs(fft(dechirp));
% plot(dechirp_fft,'-','LineWidth',6); hold on;
% plot(301,4096,'p','MarkerFaceColor','red','MarkerSize',40); hold on;
% plot([841, 1867, 3914], [460, 172, 171],'^','MarkerFaceColor','red','MarkerSize',20); hold on;
% plot([1313, 2337, 3361], [173, 571, 870],'s','MarkerFaceColor','red','MarkerSize',20); hold on;
% xlim([0 dine]);
% ylim([0 5000]);
% xlabel('BIN');
% ylabel('|FFT|'); 
% set(gca,'FontSize',80);
% set(get(gca,'YLabel'),'FontSize',100);

