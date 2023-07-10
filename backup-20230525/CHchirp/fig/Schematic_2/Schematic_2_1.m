fclose all;     %关闭所有matlab打开的文件
tic;            % 打开计时器

Config_Path = '.\config\';                                      % 设置配置文件所在路径
Verification_path = strcat(Config_Path,'sf10_BW500.txt');       % bin值验证文件
Setting_File = dir(fullfile(Config_Path,'sf10_BW500.json'));    % 配置文件
Setting_File_Path = strcat(Config_Path, Setting_File.name);
Setting_file = fopen(Setting_File_Path,'r');
setting = jsondecode(fscanf(Setting_file,'%s'));                % 解析json格式变量
lora_set.bw = setting.captures.lora_bw;                         % 设置接收数据包的lora_set.bw
lora_set.sf = setting.captures.lora_sf;                         % 设置接收数据包的lora_set.sf
lora_set.sample_rate = setting.global.core_sample_rate;         % 设置接收数据包的采样率
lora_set.sample_rate = 1e6;
lora_set.Pkg_length = setting.captures.lora_pkg_length;         % 设置接收数据包的长度
lora_set.dine = lora_set.sample_rate*bitshift(1,lora_set.sf)/lora_set.bw;    % 根据lora_set.sf和lora_set.bw计算出一个chirp包含的采样点个数
lora_set.fft_x = 2^lora_set.sf;                               % 根据lora_set.dine计算出包含lora_set.bw所需的FFT点数
lora_set.Preamble_length = 8;
leakage_width_array = [0.05,0.01,0.015,0.01];
lora_set.filter_num = 2;
lora_set.leakage_width1 = leakage_width_array(lora_set.sf-6);
lora_set.leakage_width2 = 1-lora_set.leakage_width1;

% 生成理想仿真信号 sf10,250KHz
File_Path = 'D:\CHchirp_samples\align_windows\red_node_all\samples_1.sigmf-data';
fid_1=fopen(File_Path,'rb');
[A_1]=fread(fid_1,'float32')';
A_1_length = size(A_1,2);
G0 = A_1(1:2:A_1_length-1) + A_1(2:2:A_1_length)*1i;
% 计算cfo和time_offset
[d_downchirp, d_upchirp] = build_idealchirp(lora_set, lora_set.bw/2); % build idealchirp
[cfo, windows_offset] = get_cfo_winoff(G0, lora_set, d_downchirp, d_upchirp, 2);
G0 = circshift(G0, -round(windows_offset));
[d_downchirp_cfo, d_upchirp_cfo] = rebuild_idealchirp_cfo(lora_set, cfo);
% dechirp = [d_downchirp_cfo, d_downchirp_cfo] .* G0(12.25*lora_set.dine+1:14.25*lora_set.dine);
% 画图
stft(G0(14.25*lora_set.dine+1:16.25*lora_set.dine), lora_set.sample_rate, 'Window',rectwin(64),'OverlapLength',63,'FFTLength',lora_set.dine);
ylim([-250 250]);
% imagesc(abs(s));
x = [-250 -200 -150 -100 -50 0 50 100 150 200 250];
set(gca,'YTickLabel',x);
% set(gca,'YTickLabel',str2double(get(gca,'YTickLabel'))*1000);
xlabel('Times(ms)');
ylabel('Frequency(KHz)'); 
colormap(CustomColormap);
colorbar('off');
set(gca,'FontSize',88);
set(get(gca,'XLabel'),'FontSize',88)
set(get(gca,'YLabel'),'FontSize',88)
set(gca,'YTickLabel',[-250 -125 0 125 250]);
set(gca,'xtick',[],'xticklabel',[]);
% set(gca,'xtick',[],'xticklabel',[]);
% imagesc(t2+timestart-wlen/fs/2,f,abs(tfr2))