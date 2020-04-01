void start_animation(int animation){
    if(animation > 1){
        if(animation == 3){
            std::cout << "set terminal gif animate\n"; 
            std::cout << "set output 'movie.gif'";
        }
        std::cout << "unset key\n";
        std::cout << "set xrange[-10:" << Lx+10 << "]\n";
        std::cout << "set yrange[-10:" << Ly+10 << "]\n";
        std::cout << "set size ratio -1\n";
        std::cout << "set parametric\n";
        std::cout << "set trange [0:7]\n";
        std::cout << "set isosamples 12" << std::endl;
    }
}
void begin_frame(int animation){
    if(animation > 1){
        std::cout << "plot 0,0 ";
        std::cout << " , " << Lx/7 << "*t,0";             // down wall
        std::cout << " , " << Lx/7 << "*t," << Ly;        // up wall
        std::cout << " , 0," << Ly/7 << "*t";             // left wall
        std::cout << " , " << Lx << "," << Ly/7 << "*t";  // right wall
    }
}
void end_frame(int animation){
    if (animation > 1)
        std::cout << std::endl;
}
