function SimExp_CellPlacementModel_Comparison_Plots(pred_tbl,resp_tbl,...
    str_lbl1, f1, str_lbl2, f2, pred_label,out_path)

% keyboard
model_diff=resp_tbl.(f1)-resp_tbl.(f2);

%% MCMRP to BMRP Model Comparison: Bland altman
fprintf(['Comparison of ' str_lbl1 ' to ' str_lbl2 '\n']);
% Bland Altman plot
figure;
[h,p] = ttest(resp_tbl.(f1), resp_tbl.(f2));
x=(resp_tbl.(f1)+resp_tbl.(f2))./2;
y=resp_tbl.(f1)-resp_tbl.(f2);
fprintf(['Paired T-Test ' str_lbl1 ' vs ' str_lbl2 ': %.4e\n'],p);
plot(x,y ,'ro', 'MarkerSize', 3);
hold on;
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y) mean(y)+1.96*std(y)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b'); hold off
xlabel(['\mu ICF of ' str_lbl1 ' and  ' str_lbl2]);
ylabel([str_lbl1 ' ICF - ' str_lbl2 ' ICF']);
beautifyAxis(gcf);
set(gcf,'position',[300 300 300 200])

%% MCMRP to BMRP Model Comparison: MVLR
mdl = fitlm(zscore_col(table2array(pred_tbl(:,2:end))),zscore_col(model_diff),...
    'Linear','RobustOpts','on',...
    'Intercept',true,'PredictorVars',pred_label);
mdl.Coefficients
ModelDif_Annova = anova(mdl)


%% MCMRP to BMRP Model Comparison: Pearson
fe = pred_tbl.Properties.VariableNames;fe(1)=[];
close all
clear R P RLO RUP
% Plot each predictor versus the difference between models, get R2 and CI
fprintf(['Pearson Test of [' str_lbl1 '-' str_lbl2 '] vs Predictors\n'])
for n=1:numel(fe)
    figure;
    plot(pred_tbl.(fe{n}), model_diff,'r.','MarkerSize',4)
    hold on; xa=xlim;
    plot([0 xa(2)],[0 0],'k');
    hold off;
       
    xlabel(pred_label{n});
    ylabel([str_lbl1 ' ICF - ' str_lbl2 ' ICF'])
    [R,P,RLO,RUP]= corrcoef(pred_tbl.(fe{n}), model_diff, 'alpha', 0.05);
    fprintf('%s (r=%.2f [%.2f %.2f] %s)\n',pred_label{n},...
        R(1,2),RLO(1,2),RUP(1,2),printf_pval(P(1,2)));
 
    beautifyAxis(gca); 
    set(gcf,'Position', [100 100 165 175])
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
    set(findall(gcf,'-property','FontName'),'FontName','Helvetica')
%     set(gca,'YTickLabelRotation',90)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',7)
    
    saveas(gcf,[out_path '/' str_lbl1 'vs' str_lbl2 '_' regexprep(pred_label{n},'\/','_') '.fig'])
     
    pause();
    close(gcf)
end


