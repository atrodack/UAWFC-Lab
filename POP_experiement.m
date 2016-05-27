clear all;
clc;
% close all;


voxel_num = 10^2;
POP_per_voxel = 6^2;
POP_num = voxel_num * POP_per_voxel;

FPM = zeros(sqrt(POP_num),sqrt(POP_num));
ditherer = ones(2,2);

x = linspace(1,sqrt(POP_num),sqrt(POP_num));
y = x;
[X,Y] = meshgrid(x,y);

ditherer_x = linspace(1,sqrt(POP_per_voxel),sqrt(POP_per_voxel));
ditherer_y = ditherer_x;
[ditherer_X, ditherer_Y] = meshgrid(ditherer_x,ditherer_y);


%% Plot Initial
figure(1);
subplot(1,2,1)
imagesc(x,y,FPM);
hold on
plot(X,Y,'k*');
hold off
axis xy;
axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
sqar;
subplot(1,2,2)
imagesc(ditherer_x,ditherer_y,ditherer)
hold on
plot(ditherer_X,ditherer_Y,'k*')
hold off
axis xy;
sqar;
colormap(gray(256));
drawnow;
% 
% 
% %% Set the first column
% 
% for m = 1:length(x)-1
%     for n = 1:sqrt(POP_per_voxel)
%         FPM(m,n) = FPM(m,n) + 1;
%         FPM(m+1,n) = FPM(m+1,n) + 1;
%     end
% end
% 
% figure(1);
% subplot(1,2,1)
% % imagesc(x,y,FPM);
% hold on
% for m = 1:length(x)
%     for n = 1:length(y)
%         if FPM(m,n) == 1
%             plot(X(m,n),Y(m,n),'g*');
%         elseif FPM(m,n) == 2
%             plot(X(m,n),Y(m,n),'b*');
%         elseif FPM(m,n) == 0
%             plot(X(m,n),Y(m,n),'k*');
%         end
%     end
% end
% 
% hold off
% axis xy;
% axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
% sqar;
% drawnow;
% %% Move to second column
% 
% for m = 1:length(x)-1
%     for n = 1:sqrt(POP_per_voxel)
%         FPM(m,n+1) = FPM(m,n+1) + 1;
%         FPM(m+1,n+1) = FPM(m+1,n+1) + 1;
%     end
% end
% 
% figure(1);
% subplot(1,2,1)
% % imagesc(x,y,FPM);
% hold on
% for m = 1:length(x)
%     for n = 1:length(y)
%         if FPM(m,n) == 1
%             plot(X(m,n),Y(m,n),'g*');
%         elseif FPM(m,n) == 2
%             plot(X(m,n),Y(m,n),'b*');
%         elseif FPM(m,n) == 4
%             plot(X(m,n),Y(m,n),'r*');
%         elseif FPM(m,n) == 6
%             plot(X(m,n),Y(m,n),'y*');
%         elseif FPM(m,n) == 0
%             plot(X(m,n),Y(m,n),'k*');
%         end
%     end
% end
% 
% hold off
% axis xy;
% axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
% sqar;
% drawnow;
% %% Move to third column
% for m = 1:length(x)-1
%     for n = 1:sqrt(POP_per_voxel)
%         FPM(m,n+2) = FPM(m,n+2) + 1;
%         FPM(m+1,n+2) = FPM(m+1,n+2) + 1;
%     end
% end
% 
% 
% figure(1);
% subplot(1,2,1)
% % imagesc(x,y,FPM);
% hold on
% for m = 1:length(x)
%     for n = 1:length(y)
%         if FPM(m,n) == 1
%             plot(X(m,n),Y(m,n),'g*');
%         elseif FPM(m,n) == 2
%             plot(X(m,n),Y(m,n),'b*');
%         elseif FPM(m,n) == 4
%             plot(X(m,n),Y(m,n),'r*');
%         elseif FPM(m,n) == 6
%             plot(X(m,n),Y(m,n),'y*');
%         elseif FPM(m,n) == 8
%             plot(X(m,n),Y(m,n),'c*');
%         elseif FPM(m,n) == 10
%             plot(X(m,n),Y(m,n),'m*');
%         elseif FPM(m,n) == 0
%             plot(X(m,n),Y(m,n),'k*');
%         end
%     end
% end
% 
% hold off
% axis xy;
% axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
% sqar;
% drawnow;
% 
% 
%% Full Set
clf;
FPM = zeros(sqrt(POP_num),sqrt(POP_num));
figure(1);
subplot(1,2,1)
imagesc(x,y,FPM);
hold on
plot(X,Y,'k*');
hold off
axis xy;
axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
sqar;
subplot(1,2,2)
imagesc(ditherer_x,ditherer_y,ditherer)
hold on
plot(ditherer_X,ditherer_Y,'k*')
hold off
axis xy;
sqar;
colormap(gray(256));
drawnow;


for k = 0:length(x)-sqrt(POP_per_voxel)
    for m = 1:length(x)-1
        for n = 1:sqrt(POP_per_voxel)
            FPM(m,n+k) = FPM(m,n+k) + 1;
            FPM(m+1,n+k) = FPM(m+1,n+k) + 1;             
        end
    end  
figure(1);
subplot(1,2,1)
% imagesc(x,y,FPM);
hold on
for m = 1:length(x)
    for n = k+1:k+sqrt(POP_per_voxel)
        if FPM(m,n) == 1
            plot(X(m,n),Y(m,n),'g*');
        elseif FPM(m,n) == 2
            plot(X(m,n),Y(m,n),'b*');
        elseif FPM(m,n) == 3
            plot(X(m,n),Y(m,n),'r.');
        elseif FPM(m,n) == 4
            plot(X(m,n),Y(m,n),'r*');
        elseif FPM(m,n) == 5
            plot(X(m,n),Y(m,n),'y.');
        elseif FPM(m,n) == 6
            plot(X(m,n),Y(m,n),'y*');
        elseif FPM(m,n) == 7
            plot(X(m,n),Y(m,n),'c.');
        elseif FPM(m,n) == 8
            plot(X(m,n),Y(m,n),'c*');
        elseif FPM(m,n) == 9
            plot(X(m,n),Y(m,n),'m.');
        elseif FPM(m,n) == 10
            plot(X(m,n),Y(m,n),'m*');
        elseif FPM(m,n) == 11
            plot(X(m,n),Y(m,n),'k.');
        elseif FPM(m,n) == 12
            plot(X(m,n),Y(m,n),'w*');
        elseif FPM(m,n) == 0
            plot(X(m,n),Y(m,n),'k*');
        end
    end
end

hold off
axis xy;
axis([0,sqrt(POP_num)+1,0,sqrt(POP_num)+1]);
sqar;
drawnow;
end
