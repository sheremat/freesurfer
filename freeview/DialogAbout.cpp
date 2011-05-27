/**
 * @file  DialogAbout.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/16 22:07:51 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "DialogAbout.h"
#include "ui_DialogAbout.h"

DialogAbout::DialogAbout(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogAbout)
{
  ui->setupUi(this);
  QString strg = ui->labelVersion->text();
  strg.replace("xxx", QString("%1 %2").arg(__DATE__).arg(__TIME__));
  ui->labelVersion->setText(strg);
}

DialogAbout::~DialogAbout()
{
  delete ui;
}
