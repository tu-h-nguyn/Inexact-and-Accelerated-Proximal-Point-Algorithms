function export_fig_pdf(fig, filename)
% EXPORT_FIG_PDF Ghi hinh ra PDF vector, chay duoc tren ca MATLAB lan Octave.
%
% exportgraphics() chi co tu MATLAB R2020a va khong co trong Octave, nen phai
% co duong du phong; neu khong, script chi chay duoc tren dung mot moi truong
% va CI khong the ve lai hinh.
if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
    exportgraphics(fig, filename, 'ContentType', 'vector');
else
    % Octave: print() voi thiet bi pdf cho ket qua vector tuong duong.
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, filename, '-dpdf', '-color');
end
end
