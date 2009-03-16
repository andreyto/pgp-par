#!/usr/bin/env python
# Wrapper for the libsvm library; modified slightly from the "easy.py" in the distribtion.
import sys
import os

def Main(train_pathname, test_pathname = None):
    # svm, grid, and gnuplot executable
    is_win32 = (sys.platform == 'win32')
    if not is_win32:
        svmscale_exe = "../svm-scale"
        svmtrain_exe = "../svm-train"
        svmpredict_exe = "../svm-predict"
        grid_py = "./grid.py"
        gnuplot_exe = "/usr/bin/gnuplot"
    else:
        # example for windows
        RootDir = r"C:\libsvm"
        svmscale_exe = os.path.join(RootDir,"windows","svmscale.exe")
        svmtrain_exe = os.path.join(RootDir,"windows","svmtrain.exe")
        svmpredict_exe = os.path.join(RootDir,"windows","svmpredict.exe")
        gnuplot_exe = r"c:\gnuplot\pgnuplot.exe"
        grid_py = os.path.join(RootDir,r"C:\libsvm\tools\GridNoQ.py")

    assert os.path.exists(svmscale_exe),"svm-scale executable not found"
    assert os.path.exists(svmtrain_exe),"svm-train executable not found"
    assert os.path.exists(svmpredict_exe),"svm-predict executable not found"
    #assert os.path.exists(gnuplot_exe),"gnuplot executable not found"
    assert os.path.exists(grid_py),"grid.py not found"

    assert os.path.exists(train_pathname),"training file not found"
    file_name = os.path.split(train_pathname)[1]
    scaled_file = file_name + ".scale"
    model_file = file_name + ".model"
    range_file = file_name + ".range"

    if test_pathname:
        file_name = os.path.split(test_pathname)[1]
        assert os.path.exists(test_pathname),"testing file not found"
        scaled_test_file = file_name + ".scale"
        predict_test_file = file_name + ".predict"

    cmd = "%s -s %s %s > %s" % (svmscale_exe, range_file, train_pathname, scaled_file)
    print 'Scaling training data...'
    os.system(cmd)

    cmd = "%s -svmtrain %s -gnuplot %s %s" % (grid_py, svmtrain_exe, gnuplot_exe, scaled_file)
    print 'Cross validation...'
    print "(The command is '%s')"%cmd
    dummy, f = os.popen2(cmd)

    line = ''
    while 1:
        last_line = line
        line = f.readline()
        if not line: break
    c,g,rate = map(float,last_line.split())

    print 'Best c=%s, g=%s CV rate=%s' % (c,g,rate)

    cmd = "%s -c %s -g %s %s %s" % (svmtrain_exe,c,g,scaled_file,model_file)
    print 'Training...'
    os.popen(cmd)

    print 'Output model: %s' % model_file
    if test_pathname:
        cmd = "%s -r %s %s > %s" % (svmscale_exe, range_file, test_pathname, scaled_test_file)
        print 'Scaling testing data...'
        os.system(cmd)
        
        cmd = "%s %s %s %s" % (svmpredict_exe, scaled_test_file, model_file, predict_test_file)
        print "Predict: '%s'"%cmd
        print 'Testing...'
        os.system(cmd)

        print 'Output prediction: %s' % predict_test_file

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print 'Usage: %s training_file [testing_file]' % sys.argv[0]
        raise SystemExit
    if len(sys.argv)>2:
        Main(sys.argv[1], sys.argv[2])
    else:
        Main(sys.argv[1])

