import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import functions_project as myf
import glob
import time
import function_DeltaF as dff
import os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from pandas import DataFrame as dfram
clear = lambda: os.system( 'cls' )
clear( )

class act_cells:
    '''
    1. Cell selection based on Standard deviation of the
    cells activities (observed fluorescence or df/f).
    2. Remove duplicated cells (cells with similar coordinates)
    Data structure:
         data: (signal, cell), coor: (CELL, z,y,x)

    Mahdi Zarei, 05-15-2017
    NOTES:
    >> file_coor = ('A1_merge_zyx.npz')
    >> file_fluor = ('A1_merge_fluo.npz')
    >> data= pac.act_cells(file_coor, file_fluor, no_cells=100, no_signals2show= 20, show_signals='no')
    print (np.shape(data.data))
    (500L, 8416L)
    >> print (np.min(data.dataframe.std()))
    0.0114569555514
    >> print (np.shape(data.dataframe))
    (500, 100)

    # endregion
    '''

    def __init__(self, file_coor, file_fluor, no_cells='all_cells', thr_value_std=0, calculate_cdf ='no', cell_no_steps=1, signal_start_point=0):
        coor = glob.glob(file_coor)
        fluor = glob.glob(file_fluor)

        data_input = myf.npz_coor_data2array_coor_data_array(coor, fluor)
        data = data_input.curr_data [signal_start_point: , :]
        coor = data_input.curr_coor
        num_cell_each_slice = data_input.num_cell_each_slice

        coor[:, 0] = 0.295 * coor[:, 0]
        coor[:, 1] = 0.295 * coor[:, 1]


        if no_cells == 'all_cells':
            DataFrameSize = len(data.T)
        else:
            DataFrameSize = no_cells

        df,df_coor = dff.df_cells(data,coor, num_cell= DataFrameSize, span=30, cell_no_steps=cell_no_steps)
        df = np.asarray(df)  # datafram = dfram(df.T[:, 0:])
        datafram = dfram(df.T[:, 0:DataFrameSize])

        count_thr_std = 0
        df_thr =[]
        df_thr_coor = []
        for i in range(0, len(df)):
            if datafram[i].std() > thr_value_std:
                count_thr_std += 1
                df_thr.append(np.asarray(df[i]))
                df_thr_coor.append(np.asarray(df_coor[i]))

        self.data = data_input.curr_data # Raw fluorescence
        self.coor = data_input.curr_coor
        self.num_cell_first_slice = num_cell_each_slice
        self.df=df
        self.df_coor = df_coor
        self.dataframe = datafram
        self.count_thr_std = count_thr_std
        self.df_input_size= len(df)
        self.no_cells =DataFrameSize
        self.datafram=datafram
        self.thr_value_std= thr_value_std
        self.df_thr = df_thr # OUTPUT 1
        self.df_thr_coor = df_thr_coor # OUTPUT 2

    def show_active_signal(self, cells2Show, threshold= 0, out_path= 'c:\\', file_name=' ' ):
        fig_dfc = plt.figure()
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))
        #fig_dfc.suptitle('Threshold value (STD) =' + str(threshold) + '\n\nObserved fluorescence signals             '
        #                                                              '                 DF/F' )
        try:
            os.stat(out_path + 'active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        except:
            os.mkdir(out_path + 'active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')

        ctn = 1
        for i in range(1, cells2Show):
            if self.datafram[i].std() > threshold:
                ax_ge = fig_dfc.add_subplot(cells2Show/3, 5, ctn)
             #   ax_ge.set_ylim(-0.05, 0.8)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(self.df[i])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn +=1
        plt.savefig(
            out_path + 'active_cell_'+ str(time.strftime("%Y_%m_%d_%I_%M"))+'\\' + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I")) + '.jpg')


        '''
        if calculate_cdf == 'yes':
            fig_st = plt.figure()
            fig_st.patch.set_facecolor('white')
            ax_st = fig_st.add_subplot(511)
            ax_st.plot(datafram.std(), color='red')
            plt.xlabel('CELL')
            plt.ylabel('STD')
            ax_st = fig_st.add_subplot(512)
            ax_st.plot(datafram.min())
            plt.xlabel('CELL')
            plt.ylabel('MIN')
            ax_st = fig_st.add_subplot(513)
            ax_st.plot(datafram.max())
            plt.grid('on')
            plt.xlabel('CELL')
            plt.ylabel('MAX')
            ax_st = fig_st.add_subplot(514)
            ax_st.plot(datafram.mean())
            plt.grid('on')
            plt.xlabel('TIME')
            plt.ylabel('MEAN')
            ax_st = fig_st.add_subplot(515)
            ax_st.plot(data[:, 2])
            plt.grid('on')
            plt.xlabel('Time')
            plt.ylabel('Obs. Fluo.')
            # </editor-fold>

            # <editor-fold desc="Histogram and CUM of the data">
            fig_hist = plt.figure()
            fig_hist.patch.set_facecolor('white')
            ax_hist = fig_hist.add_subplot(131)
            ax_hist.set_xlabel('STANDARD DEVIATION')
            ax_hist.set_ylabel('FREQUENCY')  # ax_hist.set_ylabel('Probability density')
            ax_hist.set_title('Histogram of fluorescence standard deviation')
            n, bins, patches = ax_hist.hist(datafram.std(), bins=len(datafram.std()))
            # values, base= np.histogram(datafram.std(), bins=40)
            # len(datafram.std())
            ax_hist_cum = fig_hist.add_subplot(132)
            # values, base = np.histogram(datafram.std(), bins=240)
            plt.plot(n, '.')
            ax_hist_cum.set_xlabel('BINS')
            ax_hist_cum.set_ylabel('FREQUENCY')  # ax_hist.set_ylabel('Probability density')
            ax_hist_cum.set_title('Histogram of fluorescence standard deviation')

            ax_cum = fig_hist.add_subplot(133)
            # cumulative_n= np.cumsum (datafram.std())
            # plt.plot(cumulative_n)
            # fig_test = plt.figure()
            # evaluate the histogram
            # values, base = np.histogram(datafram.std(), bins= 40)
            # values, base = np.histogram(datafram.std(), bins= 40)
            # evaluate the cumulative
            cumulative = np.cumsum(n)  # values
            # plot the cumulative function
            plt.plot(bins[:-1], cumulative, '.', c='green')  # plt.plot(base[:-1], cumulative, '.', c='green')
            ax_cum.set_ylabel('COUNT')
            ax_cum.set_xlabel('STANDARD DEVIATION')
            ax_cum.set_title('COMULATIVE (green)and SURVIVAL (blue) distribution')
            # plot the survival function
            plt.plot(bins[:-1], len(datafram.std()) - cumulative, '.', c='blue')
            plt.grid(True)

            sig_std = []

            for k in range(10):
                sig_std.append(np.std(df[k]))

            std_div = 10
            std_tole = (np.max(sig_std) - np.min(sig_std)) / std_div

            # </editor-fold>

            plt.show()

'''

    def show_non_active_signal(self, cells2Show, threshold=0, out_path='c:\\', file_name=' '):
        fig_dfc = plt.figure()
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))
        try:
            os.stat(out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        except:
            os.mkdir(out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        ctn = 1
        for i in range(1, cells2Show):
            if self.datafram[i].std() < threshold:
                ax_ge = fig_dfc.add_subplot(cells2Show / 3, 5, ctn)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(self.df[i])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn += 1
        plt.savefig(
            out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\' + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I")) + '.jpg')

    def show_factorAnalyis_signals(self, data, scor, cells2Show=0, threshold=0, out_path='c:\\', file_name=' '):
        '''

        :param data:
        :param scor:
        :param cells2Show:
        :param threshold:
        :param out_path:
        :param file_name:
        :return:

        Mahdi Zarei: 04-12-2017
        Example:
        >> data_input.show_factorAnalyis_signals(data_input.df, fa_scores, out_path=out_path, file_name=input_data_file)
        '''

        out_path = out_path + 'FC_cell_' + str(time.strftime("%Y_%m_%d_%I")) + '\\'
        print (out_path)
        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        #fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))

        if cells2Show == 0:
            cells2Show = len(data)

        for data_interval in range(0,len(data)/cells2Show):
            ctn = 1
            fig_dfc = plt.figure()
            fig_dfc.set_facecolor('white')
            fig_dfc.patch.set_facecolor('white')
            fig_dfc.suptitle(file_name)

            for k in range(data_interval*cells2Show, (data_interval+1)*cells2Show ):


                ax_ge = fig_dfc.add_subplot(cells2Show/3.2, 4, ctn+1)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(data[scor[k]])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn +=1

            plt.savefig(
                out_path            + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I_%M_%S")) + '.jpg', dpi=900)
            #print "~",
            plt.close()


    def imshow_factorAnalyis_signals(self, data, scor, cells2Show=100,  threshold=0, out_path='c:\\', file_name=' '):
        out_path = out_path + 'FC_cell_' + str(time.strftime("%Y_%m_%d_%I")) + '\\'
        print (out_path)
        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        # fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))


        #for data_interval in range(0, len(data) / cells2Show):
        ctn = 1
        fig_dfc = plt.figure(figsize=(7, 10))
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle(file_name)

        imgplot = plt.imshow(data[scor[1:cells2Show]],clim=(0, 0.15))
        imgplot.set_cmap('Blues_r') # nipy_spectral

        plt.colorbar(imgplot, fraction=0.048, pad=0.1)
        plt.xlabel('Time')
        plt.ylabel('CELL')


        plt.savefig(
            out_path + file_name + '_' + str(threshold) + '_' + str(
                time.strftime("%Y_%m_%d_%I_%M_%S")) + '.jpg', dpi=100)
        #print "~",
        plt.close()
#class cellSelection_factorAnalysis:
#    def __init__(self, data, coor):
#        fa=FactorAnalysis(data)












class plot_correlation_map:
    def __init__(self, corr_data, out_path='c:\\', file_name=' ', st=' ', percentage=' ', save='no', title=' '):
        fig_corr = plt.figure()
        ax1 = fig_corr.add_subplot(111)
        cmap = cm.get_cmap('jet', 30)
        cax1 = ax1.imshow(corr_data, interpolation="nearest", cmap=cmap)

        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        if save=='yes':
            #plt.savefig(
            #    out_path + file_name + st + '_' + str(percentage) + '_' + str(time.strftime("%Y_%m_%d_%I")) + '.jpg')
            #plt.savefig('Z_'+ str(time.strftime("%Y_%m_%d_%I%M"))+'.jpg', bbox_inches='tight')
            plt.title('  '+ str(title))
            plt.xticks( rotation='vertical')
            plt.savefig(
                out_path + file_name  + '.jpg', bbox_inches='tight', dpi=200)
        self.img= plt.imshow(corr_data)

            #plt.show()
def create_video( ims_video_corr, out_path , input_data_file):
    '''
    Mahdi Zarei, 06-01-2017
    Example:

        >> ims = []
        >> for fram in range(10):
                corr =...
                ims.append([plt.imshow(corr_data)])
        >> cac.create_video(ims, out_path, input_data_file)
    '''
    fig_video = plt.figure(1)
    plt.figure(1)
    ani = animation.ArtistAnimation(plt.figure(1), ims_video_corr)
    ani.save(out_path + input_data_file  + '.mp4')

class visulize_modules:
    '''
    MODULE3 THAT THEIR SIZE ARE SMALLER THAN minComSize (2) WILL NOT SHOWN
    Elements_index is the indices of the elements in each module
    '''
    colors = np.asarray(['Aqua', 'Aquamarine', 'Black', 'Blue', 'Brown', 'BurlyWood', 'CadetBlue',
                         'Chartreuse', 'Chocolate', 'Coral', 'CornflowerBlue', 'Crimson', 'Cyan',
                         'DarkBlue', 'DarkCyan', 'DarkGoldenRod', 'DarkGrey', 'DarkGreen', 'DarkKhaki',
                         'DarkMagenta', 'DarkOliveGreen', 'Darkorange', 'DarkRed', 'DarkSalmon',
                         'DarkSeaGreen', 'DarkSlateBlue', 'DarkSlateGray', 'DarkTurquoise',
                         'DarkViolet', 'DeepPink', 'DeepSkyBlue', 'DimGray', 'DimGrey', 'FireBrick',
                         'ForestGreen', 'Fuchsia', 'Gold', 'GoldenRod', 'Gray', 'Grey', 'Green',
                         'GreenYellow', 'HotPink', 'IndianRed,', 'Indigo,', 'Khaki', 'LawnGreen',
                         'LemonChiffon', 'LightBlue', 'LightCoral', 'LightGray', 'LightGrey',
                         'LightGreen', 'LightPink', 'LightSalmon', 'LightSeaGreen', 'LightSkyBlue',
                         'LightSlateGray', 'Lime', 'LimeGreen', 'Magenta', 'Maroon', 'MediumAquaMarine',
                         'MediumBlue', 'MediumOrchid', 'MediumPurple', 'MediumSeaGreen',
                         'MediumSlateBlue', 'MediumSpringGreen', 'MediumTurquoise', 'MediumVioletRed',
                         'MidnightBlue', 'MintCream', 'MistyRose', 'Moccasin', 'NavajoWhite', 'Navy',
                         'OldLace', 'Olive', 'OliveDrab', 'Orange', 'OrangeRed', 'Orchid',
                         'PaleGoldenRod', 'PaleGreen', 'PaleTurquoise', 'PaleVioletRed', 'PapayaWhip',
                         'PeachPuff', 'Peru', 'Pink', 'Plum', 'PowderBlue', 'Purple', 'Red',
                         'RosyBrown', 'RoyalBlue', 'SaddleBrown', 'Salmon', 'SandyBrown', 'SeaGreen',
                         'SeaShell', 'Sienna', 'Silver', 'SkyBlue', 'SlateBlue', 'SlateGray',
                         'SlateGrey', 'Snow', 'SpringGreen', 'SteelBlue', 'Tan', 'Teal', 'Thistle',
                         'Tomato', 'Turquoise', 'Violet', 'Wheat', 'White', 'WhiteSmoke', 'Yellow',
                         'YellowGreen'])
    CurrCommunity_index = 0
    nuclei_size =15
    step_v = 15
    size_v = 2

    def __init__(self, df_coord, df, corr, ci, minComSize = 1):

        df_coord = np.asarray(df_coord)
        ci = np.asarray(ci)
        df = np.asarray(df)
        module_coord, module_df, module_corr =[],[], []
        CurrCommunity_index = 1
        no_community = 0
        #minComSize = 1
        for module_index in range(1, max(ci)):
            elements_index = myf.find_index(np.asarray(ci), module_index)
            if len(elements_index) > minComSize:
                no_community +=1
        for module_index in range(1, max(ci)):
            elements_index = myf.find_index(np.asarray(ci), module_index)
            module_coord.append(df_coord[elements_index])
            module_df.append(df[elements_index])
            mo_corr_temp = corr[elements_index]
            module_corr.append(mo_corr_temp[:,elements_index])

        self.no_communities = no_community
        self.module_coord = module_coord
        self.module_df = module_df
        self.module_corr= module_corr
        self.ci=ci
        self.coord= df_coord
        self.df=df

    # Visualization of the modules
    def visulize_modules_together(self, coordinate, minComSize=1, no_cells=100, max_no_module2show=100):
        fig_all_module = plt.figure()
        fig_all_module.patch.set_facecolor("W")
        ax_all_module = fig_all_module.add_subplot(1, 2, 1, projection='3d')
        ax_all_module_02 = fig_all_module.add_subplot(1, 2, 2, projection='3d')
        ax_all_module_02.set_zlim([-100, 150])

        # Get rid of the panes
        ax_all_module.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Get rid of the spines
        ax_all_module.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        # Get rid of the ticks
        ax_all_module.set_xticks([])
        ax_all_module.set_yticks([])
        ax_all_module.set_zticks([])
        ax_all_module.view_init(0, -90)

        # Get rid of the panes
        ax_all_module_02.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module_02.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module_02.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Get rid of the spines
        ax_all_module_02.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module_02.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax_all_module_02.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        # Get rid of the ticks
        ax_all_module_02.set_xticks([])
        ax_all_module_02.set_yticks([])
        ax_all_module_02.set_zticks([])

        ax_all_module.scatter(coordinate[0:no_cells:self.step_v, 0],
                              coordinate[0:no_cells:self.step_v, 1],
                              coordinate[0:no_cells:self.step_v, 2], s=self.size_v,
                              marker='o', color='lightgreen', c='lightgreen')
        ax_all_module_02.scatter(coordinate[0:no_cells:self.step_v, 0],
                                 coordinate[0:no_cells:self.step_v, 1],
                                 coordinate[0:no_cells:self.step_v, 2], s=self.size_v,
                                 marker='o', color='lightgreen', c='lightgreen')

        for module_index in range(1, max(self.ci)):
            elements_index = myf.find_index(np.asarray(self.ci), module_index)
            if len(elements_index) > minComSize and module_index <max_no_module2show:
                self.CurrCommunity_index += 1
                ax_all_module.set_ylabel("  "), ax_all_module.set_xlabel("  ")  # Frontal
                ax_all_module.set_zlim([-100, 150])
                # ax_all_module.axis([100, 1050, 100, 1050])
                mngr = plt.get_current_fig_manager()
                mngr.window.setGeometry(15, 35, 1000, 500)

                for CurrentCommunity_Elements in range(0, len(elements_index)):

                    for current_community_elements in range(0, len(elements_index)):
                        module_coor_current_module = np.asarray(
                            self.coord[elements_index])  # coor[element_index[:]] is the elements coordination

                        ax_all_module.scatter(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                                              module_coor_current_module[:, 2], s=self.nuclei_size,
                                              color=self.colors[self.CurrCommunity_index],
                                              marker='o', c=self.colors[self.CurrCommunity_index])  # self.colors[module_index]

                        ax_all_module_02.scatter(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                                                 module_coor_current_module[:, 2], s=self.nuclei_size,
                                                 color=self.colors[self.CurrCommunity_index],
                                                 marker='o', c=self.colors[self.CurrCommunity_index])  # self.colors[module_index]

#       plt.show()




    def visulize_signals_and_modules (self, coordinate, minComSize=1, no_cells=100, num_cell_first_slice=300, max_no_module2show=100, file_name=' ', percentage_cells=' ', out_path='c:\\'):
        AllModules_plot = plt.figure()
        AllModules_plot.suptitle(file_name)
        AllModules_plot.patch.set_facecolor('W')
        self.no_communities =  max_no_module2show #~ 05242017
        CurrCommunity_index = 0

        for module_index in range(1, min (30, max(self.ci))):  # min (max_no_module2show, max(self.ci))
            elements_index = myf.find_index(np.asarray(self.ci), module_index)
            if len(elements_index) > minComSize and len(elements_index) < 15: #and module_index < max_no_module2show:
                CurrCommunity_index += 1

                # Show module separately
                ax_AllModules_sig = AllModules_plot.add_subplot(self.no_communities, 2, 2 * CurrCommunity_index -1)
                for CurrentCommunity_Elements in range(0, len(elements_index)):
                    ax_AllModules_sig.plot(self.df[elements_index[CurrentCommunity_Elements]])
                    ax_AllModules_sig.set_ylabel(" "), ax_AllModules_sig.set_xlabel(" TIME  ")
                    ax_AllModules_sig.set_ylim([-0.05, 1.5])
                    #ax_AllModules_sig.set_xticklabels([])
                    ax_AllModules_sig.set_yticklabels([])

                ax_Module_3d = AllModules_plot.add_subplot(self.no_communities, 6,
                                                                6 * CurrCommunity_index, projection='3d')
                ax_Module_3d.set_title('m ' + str(module_index), fontsize=7)

                ax_Module_3d.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d.view_init(0, 0)

                ax_Module_3d_v2 = AllModules_plot.add_subplot(self.no_communities, 6,
                                                                   6 * CurrCommunity_index - 1,
                                                                   projection='3d')
                ax_Module_3d_v2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v2.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v2.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v2.view_init(90, 0)

                ax_Module_3d_v3 = AllModules_plot.add_subplot(self.no_communities, 6,
                                                                   6 * CurrCommunity_index - 2,
                                                                   projection='3d')
                ax_Module_3d_v3.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v3.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v3.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v3.view_init(0, 90)

                for current_community_elements in range(0, len(elements_index)):
                    module_coor_current_module = np.asarray(
                        self.coord[elements_index])  # coor[element_index[:]] is the elements coordination

                    # ax_all_module.plot(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                    #                     module_coor_current_module[:, 2], '.')

                    ax_Module_3d.scatter(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                                         module_coor_current_module[:, 2],
                                         color=self.colors[CurrCommunity_index], s=4,
                                         marker='o', c=self.colors[CurrCommunity_index])

                    ax_Module_3d.scatter(coordinate[1:num_cell_first_slice:2, 0],
                                         coordinate[1:num_cell_first_slice:2, 1],
                                         coordinate[1:num_cell_first_slice:2, 2], s=self.size_v, marker='.',
                                         color='snow')

                    ax_Module_3d_v2.scatter(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                                            module_coor_current_module[:, 2],
                                            color=self.colors[CurrCommunity_index], s=4,
                                            marker='o', c=self.colors[CurrCommunity_index])

                    ax_Module_3d_v2.scatter(coordinate[1:num_cell_first_slice:2, 0],
                                            coordinate[1:num_cell_first_slice:2, 1],
                                            coordinate[1:num_cell_first_slice:2, 2], s=self.size_v, marker='.',
                                            color='snow')

                    ax_Module_3d_v3.scatter(module_coor_current_module[:, 0], module_coor_current_module[:, 1],
                                            module_coor_current_module[:, 2],
                                            color=self.colors[CurrCommunity_index], s=4,
                                            marker='o', c=self.colors[CurrCommunity_index])

                    ax_Module_3d_v3.scatter(coordinate[1:num_cell_first_slice:2, 0],
                                            coordinate[1:num_cell_first_slice:2, 1],
                                            coordinate[1:num_cell_first_slice:2, 2], s=self.size_v, marker='.',
                                            color='snow')

                    ax_Module_3d.set_ylabel("  ", fontsize=7)

                    ax_Module_3d.axis([25, 250, 80, 225])
                    ax_Module_3d_v2.axis([65, 250, 80, 225])
                    ax_Module_3d_v3.axis([65, 250, 80, 225])

                #plt.axis([200, 720, 150, 750])
                #ax_Module_3d.set_yticklabels([])
                #ax_Module_3d.set_xticklabels([])
                #ax_Module_3d.set_zticklabels([])
                #plt.xlim([-50, 1000])
                #plt.ylim([-100, 1100])

                plt.tight_layout()
                # ZOOM
                # ax_all_module.axis([230*0.295, 800*0.295, 300*0.295, 730*0.295])
                # Get rid of the spines
                ax_Module_3d.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                # Get rid of the ticks
                ax_Module_3d.set_xticks([])
                ax_Module_3d.set_yticks([])
                ax_Module_3d.set_zticks([])

                # Get rid of the spines
                ax_Module_3d_v2.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v2.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v2.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                # Get rid of the ticks
                ax_Module_3d_v2.set_xticks([])
                ax_Module_3d_v2.set_yticks([])
                ax_Module_3d_v2.set_zticks([])

                # Get rid of the spines
                ax_Module_3d_v3.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v3.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                ax_Module_3d_v3.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
                # Get rid of the ticks
                ax_Module_3d_v3.set_xticks([])
                ax_Module_3d_v3.set_yticks([])
                ax_Module_3d_v3.set_zticks([])

                ax_Module_3d.set_xlabel('P')
                #ax_Module_3d_v2.set_xlabel('P')
                #ax_Module_3d_v3.set_xlabel('P')

        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(900, 35, 1000, 980)
        AllModules_plot.tight_layout()
        plt.savefig(out_path+file_name + '_' + str(percentage_cells)+ '_' + str(time.strftime("%Y_%m_%d_%I"))+ '.jpg' )
        #plt.show()



class file2df_coor:
    '''
    1. Cell selection based on Standard deviation of the
    cells activities (observed fluorescence or df/f).
    2. Remove duplicated cells (cells with similar coordinates)
    Data structure:
         data: (signal, cell), coor: (CELL, z,y,x)

    Mahdi Zarei, UCSF, 05-15-2017
    NOTES:
    >> file_coor = ('A1_merge_zyx.npz')
    >> file_fluor = ('A1_merge_fluo.npz')
    >> data= pac.act_cells(file_coor, file_fluor, no_cells=100, no_signals2show= 20, show_signals='no')
    print (np.shape(data.data))
    (500L, 8416L)
    >> print (np.min(data.dataframe.std()))
    0.0114569555514
    >> print (np.shape(data.dataframe))
    (500, 100)

    # endregion
    '''

    def __init__(self, filePath , no_cells='all_cells', thr_value_std=0, calculate_cdf ='no', cell_no_steps=1, signal_start_point=0):
        #coor = glob.glob(file_coor)
        #fluor = glob.glob(file_fluor)
        #files = glob.glob (filePath)

        #data_input = myf.npz_coor_data2array_coor_data_array(coor, fluor)
        #data_input = myf.readFiles (filePath)
        data_input = myf.readFiles_merge(filePath)
        data = data_input.fluor [signal_start_point: , :]
        coor = data_input.coor
        num_cell_each_slice = data_input.num_cell_each_slice

        coor[:, [0, 2]] = coor[:, [2, 0]]

        if no_cells == 'all_cells':
            DataFrameSize = len(data.T)
        else:
            DataFrameSize = no_cells

        df,df_coor = dff.df_cells(data,coor, num_cell= DataFrameSize, span=30, cell_no_steps=cell_no_steps)
        df = np.asarray(df)  # datafram = dfram(df.T[:, 0:])
        datafram = dfram(df.T[:, 0:DataFrameSize])

        corrMat = np.abs(np.corrcoef(df))
        count_thr_std = 0
        df_thr =[]
        df_thr_coor = []
        for i in range(0, len(df)):
            if datafram[i].std() > thr_value_std:
                count_thr_std += 1
                df_thr.append(np.asarray(df[i]))
                df_thr_coor.append(np.asarray(df_coor[i]))

        self.data = data_input.fluor # Raw fluorescence
        self.coor = data_input.coor
        self.num_cell_first_slice = num_cell_each_slice
        self.df=df
        self.df_coor = df_coor
        self.dataframe = datafram
        self.count_thr_std = count_thr_std
        self.df_input_size= len(df)
        self.no_cells =DataFrameSize
        self.datafram=datafram
        self.thr_value_std= thr_value_std
        self.df_thr = df_thr # OUTPUT 1
        self.df_thr_coor = df_thr_coor # OUTPUT 2
        self.df_corrMat = np.abs(np.corrcoef(corrMat))

    def show_active_signal(self, cells2Show, threshold= 0, out_path= 'c:\\', file_name=' ' ):
        fig_dfc = plt.figure()
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))
        #fig_dfc.suptitle('Threshold value (STD) =' + str(threshold) + '\n\nObserved fluorescence signals             '
        #                                                              '                 DF/F' )
        try:
            os.stat(out_path + 'active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        except:
            os.mkdir(out_path + 'active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')

        ctn = 1
        for i in range(1, cells2Show):
            if self.datafram[i].std() > threshold:
                ax_ge = fig_dfc.add_subplot(cells2Show/3, 5, ctn)
             #   ax_ge.set_ylim(-0.05, 0.8)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(self.df[i])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn +=1
        plt.savefig(
            out_path + 'active_cell_'+ str(time.strftime("%Y_%m_%d_%I_%M"))+'\\' + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I")) + '.jpg')

    def show_non_active_signal(self, cells2Show, threshold=0, out_path='c:\\', file_name=' '):
        fig_dfc = plt.figure()
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))
        try:
            os.stat(out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        except:
            os.mkdir(out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\')
        ctn = 1
        for i in range(1, cells2Show):
            if self.datafram[i].std() < threshold:
                ax_ge = fig_dfc.add_subplot(cells2Show / 3, 5, ctn)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(self.df[i])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn += 1
        plt.savefig(
            out_path + 'non_active_cell_' + str(time.strftime("%Y_%m_%d_%I_%M")) + '\\' + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I")) + '.jpg')

    def show_factorAnalyis_signals(self, data, scor, cells2Show=0, threshold=0, out_path='c:\\', file_name=' '):
        '''

        :param data:
        :param scor:
        :param cells2Show:
        :param threshold:
        :param out_path:
        :param file_name:
        :return:

        Mahdi Zarei: 04-12-2017
        Example:
        >> data_input.show_factorAnalyis_signals(data_input.df, fa_scores, out_path=out_path, file_name=input_data_file)
        '''

        out_path = out_path + 'FC_cell_' + str(time.strftime("%Y_%m_%d_%I")) + '\\'
        print (out_path)
        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        #fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))

        if cells2Show == 0:
            cells2Show = len(data)

        for data_interval in range(0,len(data)/cells2Show):
            ctn = 1
            fig_dfc = plt.figure()
            fig_dfc.set_facecolor('white')
            fig_dfc.patch.set_facecolor('white')
            fig_dfc.suptitle(file_name)

            for k in range(data_interval*cells2Show, (data_interval+1)*cells2Show ):


                ax_ge = fig_dfc.add_subplot(cells2Show/3.2, 4, ctn+1)
                ax_ge.set_ylim(-0.05, 0.17)
                cax_ge = plt.plot(data[scor[k]])
                plt.xlabel('Time')
                plt.ylabel('OBS. FLU.')
                plt.axis('off')
                ctn +=1

            plt.savefig(
                out_path            + file_name + '_' + str(threshold) + '_' + str(time.strftime("%Y_%m_%d_%I_%M_%S")) + '.jpg', dpi=900)
            #print "~",
            plt.close()


    def imshow_factorAnalyis_signals(self, data, scor, cells2Show=100,  threshold=0, out_path='c:\\', file_name=' '):
        out_path = out_path + 'FC_cell_' + str(time.strftime("%Y_%m_%d_%I")) + '\\'
        print (out_path)
        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        # fig_dfc.suptitle('Threshold value (STD) =' + str(threshold))


        #for data_interval in range(0, len(data) / cells2Show):
        ctn = 1
        fig_dfc = plt.figure(figsize=(7, 10))
        fig_dfc.set_facecolor('white')
        fig_dfc.patch.set_facecolor('white')
        fig_dfc.suptitle(file_name)

        imgplot = plt.imshow(data[scor[1:cells2Show]],clim=(0, 0.15))
        imgplot.set_cmap('Blues_r') # nipy_spectral

        plt.colorbar(imgplot, fraction=0.048, pad=0.1)
        plt.xlabel('Time')
        plt.ylabel('CELL')
        fc_fig_dir = out_path + file_name + '_' + str(threshold) + '_' + str(
            time.strftime("%Y_%m_%d_%I_%M_%S")) + '.jpg'
        plt.savefig(fc_fig_dir, dpi=100)
        plt.close()
        self.fc_fig_dir = fc_fig_dir

    