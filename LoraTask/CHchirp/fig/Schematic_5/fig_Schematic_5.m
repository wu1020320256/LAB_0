tmp_result = zeros(4,10);
for sf = 7:10
    file_name = 'result_array_sf' + string(sf) + '.mat';
    exampleObject = matfile(file_name);
    varlist = who(exampleObject);
    for index = 2:2
        varName = varlist{index};
        result_array = exampleObject.(varName);
        tmp_result(sf-6,6:10) = result_array(1, 1:2:9)./result_array(2, 1:2:9);
        tmp_result(sf-6,1:5) = result_array(1, 10:-2:2)./result_array(2, 10:-2:2);
    end
end
% tmp_result = zeros(1,10);
% tmp_result(6:10) = result_array(1, 1:2:9)./result_array(2, 1:2:9);
% tmp_result(1:5) = result_array(1, 10:-2:2)./result_array(2, 10:-2:2);
tmptmp_result = tmp_result';
bar(tmptmp_result);
set(gca,'XTickLabel', {'-434.175','-433.975','-433.775','-433.575','-433.375','433.375','433.575','433.775','433.975','434.175'});
legend('sf7','sf8','sf9','sf10');