function plot_prob(k,prev,dist,color)

%         title(sprintf('time = %d',k));
%         subplot(311); 
%         histfit(a_rand,20,'kernel'); title('a');
%         plot(x_axis,pdf(a_dist,x_axis));
%         
        plot([k-1, k],[prev,mean(dist)],color,'linewidth',2); hold on
        scatter([k, k],[max(dist.mean-3*dist.sigma,0),dist.mean+3*dist.sigma],color,'.','MarkerEdgeAlpha',.2);
        xlim([0 605]); 
%         yt = get(gca, 'YTick'); set(gca, 'YTick', yt, 'YTickLabel', yt/numel(a_rand))
        
%         subplot(312); 
%         histfit(b_rand,20,'kernel'); title('b');
%         xlim([5 25]); 
%         yt = get(gca, 'YTick'); set(gca, 'YTick', yt, 'YTickLabel', yt/numel(b_rand))
%         
%         subplot(313); histfit(hm_rand,20,'kernel'); title('hm');
%         xlim([6 25]); 
%         yt = get(gca, 'YTick'); set(gca, 'YTick', yt, 'YTickLabel', yt/numel(hm_rand))
        drawnow;
end

