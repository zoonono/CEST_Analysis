%Calculating the mean and std of the baselines 
A = zeros(1,input.number_of_baseline_sets); 
C = zeros(1,input.number_of_baseline_sets); 
for i = 1: input.number_of_baseline_sets
    A(i) = auc_STAsym_sum_roi(i).auc_around_mid_point;
    C(i) = auc_BS_sum_roi(i).auc_around_mid_point;
end
mean_auc_glucose_STAsym_sum_roi = mean(A);
 std_auc_glucose_STAsym_sum_roi = std(A);
mean_auc_glucose_BS_sum_roi     = mean(C);
 std_auc_glucose_BS_sum_roi     = std(C);

%Calculating the mean and std of the post injection scans 
A1 = zeros(1,input.number_of_post_inj_sets); 
C1 = zeros(1,input.number_of_post_inj_sets); 
for i = 1: input.number_of_post_inj_sets
    A1(i) = auc_STAsym_sum_roi(i+input.number_of_baseline_sets).auc_around_mid_point;
    C1(i) = auc_BS_sum_roi(i+input.number_of_baseline_sets).auc_around_mid_point;
end
mean_auc_glucose2_STAsym_sum_roi = mean(A1);
 std_auc_glucose2_STAsym_sum_roi = std(A1);
mean_auc_glucose2_BS_sum_roi     = mean(C1);
 std_auc_glucose2_BS_sum_roi     = std(C1);
