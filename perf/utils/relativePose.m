function rel = relativePose(curr_pose, prev_pose)

offseted = curr_pose;
offseted_prev = prev_pose;
offseted(1:3, 4) = offseted(1:3, 4) - offseted_prev(1:3, 4);
offseted_prev(1:3, 4) = 0;
rel = offseted_prev \ offseted;

end