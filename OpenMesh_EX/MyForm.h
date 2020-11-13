#pragma once

#include "DotNetUtilities.h"
#include "Mesh/GUA_OM.h"
#include "Mesh/DP.h"

Tri_Mesh *mesh;

xform xf;
GLCamera camera;
float fov = 0.7f;

static const Mouse::button physical_to_logical_map[] = {
	Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
	Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
};
Mouse::button Mouse_State = Mouse::ROTATE;

namespace OpenMesh_EX {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	bool closeFace = false;
	/// <summary>
	/// MyForm 的摘要
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		int paintStatus = 0;
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO:  在此加入建構函式程式碼
			//
		}

	protected:
		/// <summary>
		/// 清除任何使用中的資源。
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}

	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  fileToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  loadModelToolStripMenuItem;
	private: System::Windows::Forms::OpenFileDialog^  openModelDialog;
	private: System::Windows::Forms::SaveFileDialog^  saveModelDialog;
	private: System::Windows::Forms::ToolStripMenuItem^  saveModelToolStripMenuItem;
	private: HKOGLPanel::HKOGLPanelControl^  hkoglPanelControl1;


	private: System::Windows::Forms::HScrollBar^ hScrollBar2;
	private: System::Windows::Forms::Label^ Normal_Error;



	private: System::Windows::Forms::ToolStripMenuItem^ modeToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ connectColorToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ closeFaceToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ nonCloseFaceToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ faceClusterToolStripMenuItem;
	private: System::Windows::Forms::GroupBox^ Setting;
	private: System::Windows::Forms::Label^ Face_Error;

	private: System::Windows::Forms::HScrollBar^ hScrollBar1;
	private: System::Windows::Forms::Label^ label3;
	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::ToolStripMenuItem^ faceClusterCloseFaceToolStripMenuItem;

	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::NumericUpDown^ ScaleNum;
	private: System::Windows::Forms::Label^ Scale_Num;
	private: System::Windows::Forms::Label^ Ignore_Num;

	private: System::Windows::Forms::HScrollBar^ hScrollBar3;
	private: System::Windows::Forms::Label^ face_Ignore;
	private: System::Windows::Forms::ToolStripMenuItem^ saveBoundaryToolStripMenuItem;
	private: System::ComponentModel::BackgroundWorker^ backgroundWorker1;
	private: System::Windows::Forms::ToolStripMenuItem^ paintToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ originToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ newToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ saveNewModelToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ saveNewModelMtlToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ faceClusterpopulosityToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ hypothesismeshToolStripMenuItem;



	protected:

	private:
		/// <summary>
		/// 設計工具所需的變數。
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
		/// 這個方法的內容。
		/// </summary>
		void InitializeComponent(void)
		{
			HKOGLPanel::HKCOGLPanelCameraSetting^ hkcoglPanelCameraSetting1 = (gcnew HKOGLPanel::HKCOGLPanelCameraSetting());
			HKOGLPanel::HKCOGLPanelPixelFormat^ hkcoglPanelPixelFormat1 = (gcnew HKOGLPanel::HKCOGLPanelPixelFormat());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->fileToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->loadModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveBoundaryToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveNewModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveNewModelMtlToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->modeToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->connectColorToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->closeFaceToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->nonCloseFaceToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->faceClusterToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->faceClusterCloseFaceToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->paintToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->originToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->newToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openModelDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->saveModelDialog = (gcnew System::Windows::Forms::SaveFileDialog());
			this->hkoglPanelControl1 = (gcnew HKOGLPanel::HKOGLPanelControl());
			this->hScrollBar2 = (gcnew System::Windows::Forms::HScrollBar());
			this->Normal_Error = (gcnew System::Windows::Forms::Label());
			this->Setting = (gcnew System::Windows::Forms::GroupBox());
			this->Ignore_Num = (gcnew System::Windows::Forms::Label());
			this->hScrollBar3 = (gcnew System::Windows::Forms::HScrollBar());
			this->face_Ignore = (gcnew System::Windows::Forms::Label());
			this->Scale_Num = (gcnew System::Windows::Forms::Label());
			this->ScaleNum = (gcnew System::Windows::Forms::NumericUpDown());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->Face_Error = (gcnew System::Windows::Forms::Label());
			this->hScrollBar1 = (gcnew System::Windows::Forms::HScrollBar());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->backgroundWorker1 = (gcnew System::ComponentModel::BackgroundWorker());
			this->faceClusterpopulosityToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->hypothesismeshToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->menuStrip1->SuspendLayout();
			this->Setting->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ScaleNum))->BeginInit();
			this->SuspendLayout();
			// 
			// menuStrip1
			// 
			this->menuStrip1->ImageScalingSize = System::Drawing::Size(20, 20);
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(3) {
				this->fileToolStripMenuItem,
					this->modeToolStripMenuItem, this->paintToolStripMenuItem
			});
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Padding = System::Windows::Forms::Padding(4, 2, 0, 2);
			this->menuStrip1->Size = System::Drawing::Size(613, 24);
			this->menuStrip1->TabIndex = 1;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// fileToolStripMenuItem
			// 
			this->fileToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(5) {
				this->loadModelToolStripMenuItem,
					this->saveModelToolStripMenuItem, this->saveBoundaryToolStripMenuItem, this->saveNewModelToolStripMenuItem, this->saveNewModelMtlToolStripMenuItem
			});
			this->fileToolStripMenuItem->Name = L"fileToolStripMenuItem";
			this->fileToolStripMenuItem->Size = System::Drawing::Size(38, 20);
			this->fileToolStripMenuItem->Text = L"File";
			// 
			// loadModelToolStripMenuItem
			// 
			this->loadModelToolStripMenuItem->Name = L"loadModelToolStripMenuItem";
			this->loadModelToolStripMenuItem->Size = System::Drawing::Size(187, 22);
			this->loadModelToolStripMenuItem->Text = L"Load Model";
			this->loadModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::loadModelToolStripMenuItem_Click);
			// 
			// saveModelToolStripMenuItem
			// 
			this->saveModelToolStripMenuItem->Name = L"saveModelToolStripMenuItem";
			this->saveModelToolStripMenuItem->Size = System::Drawing::Size(187, 22);
			this->saveModelToolStripMenuItem->Text = L"Save Model";
			this->saveModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveModelToolStripMenuItem_Click);
			// 
			// saveBoundaryToolStripMenuItem
			// 
			this->saveBoundaryToolStripMenuItem->Name = L"saveBoundaryToolStripMenuItem";
			this->saveBoundaryToolStripMenuItem->Size = System::Drawing::Size(187, 22);
			this->saveBoundaryToolStripMenuItem->Text = L"Save Boundary";
			this->saveBoundaryToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveBoundaryToolStripMenuItem_Click);
			// 
			// saveNewModelToolStripMenuItem
			// 
			this->saveNewModelToolStripMenuItem->Name = L"saveNewModelToolStripMenuItem";
			this->saveNewModelToolStripMenuItem->Size = System::Drawing::Size(187, 22);
			this->saveNewModelToolStripMenuItem->Text = L"Save NewModel";
			this->saveNewModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveNewModelToolStripMenuItem_Click);
			// 
			// saveNewModelMtlToolStripMenuItem
			// 
			this->saveNewModelMtlToolStripMenuItem->Name = L"saveNewModelMtlToolStripMenuItem";
			this->saveNewModelMtlToolStripMenuItem->Size = System::Drawing::Size(187, 22);
			this->saveNewModelMtlToolStripMenuItem->Text = L"Save NewModelMtl";
			this->saveNewModelMtlToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveNewModelMtlToolStripMenuItem_Click);
			// 
			// modeToolStripMenuItem
			// 
			this->modeToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(7) {
				this->connectColorToolStripMenuItem,
					this->closeFaceToolStripMenuItem, this->nonCloseFaceToolStripMenuItem, this->faceClusterToolStripMenuItem, this->faceClusterCloseFaceToolStripMenuItem,
					this->faceClusterpopulosityToolStripMenuItem , this->hypothesismeshToolStripMenuItem
			});
			this->modeToolStripMenuItem->Name = L"modeToolStripMenuItem";
			this->modeToolStripMenuItem->Size = System::Drawing::Size(54, 20);
			this->modeToolStripMenuItem->Text = L"Mode";
			// 
			// connectColorToolStripMenuItem
			// 
			this->connectColorToolStripMenuItem->Name = L"connectColorToolStripMenuItem";
			this->connectColorToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->connectColorToolStripMenuItem->Text = L"Connect Color";
			this->connectColorToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::connectColorToolStripMenuItem_Click);
			// 
			// closeFaceToolStripMenuItem
			// 
			this->closeFaceToolStripMenuItem->Name = L"closeFaceToolStripMenuItem";
			this->closeFaceToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->closeFaceToolStripMenuItem->Text = L"CloseFace";
			this->closeFaceToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::closeFaceToolStripMenuItem_Click);
			// 
			// nonCloseFaceToolStripMenuItem
			// 
			this->nonCloseFaceToolStripMenuItem->Name = L"nonCloseFaceToolStripMenuItem";
			this->nonCloseFaceToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->nonCloseFaceToolStripMenuItem->Text = L"Non-CloseFace";
			this->nonCloseFaceToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::nonCloseFaceToolStripMenuItem_Click);
			// 
			// faceClusterToolStripMenuItem
			// 
			this->faceClusterToolStripMenuItem->Name = L"faceClusterToolStripMenuItem";
			this->faceClusterToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->faceClusterToolStripMenuItem->Text = L"faceCluster";
			this->faceClusterToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::faceClusterToolStripMenuItem_Click);
			// 
			// faceClusterCloseFaceToolStripMenuItem
			// 
			this->faceClusterCloseFaceToolStripMenuItem->Name = L"faceClusterCloseFaceToolStripMenuItem";
			this->faceClusterCloseFaceToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->faceClusterCloseFaceToolStripMenuItem->Text = L"faceCluster(CloseFace)";
			this->faceClusterCloseFaceToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::faceClusterCloseFaceToolStripMenuItem_Click);
			// 
			// paintToolStripMenuItem
			// 
			this->paintToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->originToolStripMenuItem,
					this->newToolStripMenuItem
			});
			this->paintToolStripMenuItem->Name = L"paintToolStripMenuItem";
			this->paintToolStripMenuItem->Size = System::Drawing::Size(47, 20);
			this->paintToolStripMenuItem->Text = L"Paint";
			// 
			// originToolStripMenuItem
			// 
			this->originToolStripMenuItem->Name = L"originToolStripMenuItem";
			this->originToolStripMenuItem->Size = System::Drawing::Size(180, 22);
			this->originToolStripMenuItem->Text = L"Origin";
			this->originToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::originToolStripMenuItem_Click);
			// 
			// newToolStripMenuItem
			// 
			this->newToolStripMenuItem->Name = L"newToolStripMenuItem";
			this->newToolStripMenuItem->Size = System::Drawing::Size(180, 22);
			this->newToolStripMenuItem->Text = L"New";
			this->newToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::newToolStripMenuItem_Click);
			// 
			// openModelDialog
			// 
			this->openModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::openModelDialog_FileOk);
			// 
			// saveModelDialog
			// 
			this->saveModelDialog->DefaultExt = L"obj";
			this->saveModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::saveModelDialog_FileOk);
			// 
			// hkoglPanelControl1
			// 
			hkcoglPanelCameraSetting1->Far = 1000;
			hkcoglPanelCameraSetting1->Fov = 45;
			hkcoglPanelCameraSetting1->Near = -1000;
			hkcoglPanelCameraSetting1->Type = HKOGLPanel::HKCOGLPanelCameraSetting::CAMERATYPE::ORTHOGRAPHIC;
			this->hkoglPanelControl1->Camera_Setting = hkcoglPanelCameraSetting1;
			this->hkoglPanelControl1->Dock = System::Windows::Forms::DockStyle::Fill;
			this->hkoglPanelControl1->Location = System::Drawing::Point(0, 24);
			this->hkoglPanelControl1->Name = L"hkoglPanelControl1";
			hkcoglPanelPixelFormat1->Accumu_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Alpha_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Stencil_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			this->hkoglPanelControl1->Pixel_Format = hkcoglPanelPixelFormat1;
			this->hkoglPanelControl1->Size = System::Drawing::Size(613, 430);
			this->hkoglPanelControl1->TabIndex = 2;
			this->hkoglPanelControl1->Load += gcnew System::EventHandler(this, &MyForm::hkoglPanelControl1_Load);
			this->hkoglPanelControl1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MyForm::hkoglPanelControl1_Paint);
			this->hkoglPanelControl1->KeyDown += gcnew System::Windows::Forms::KeyEventHandler(this, &MyForm::hkoglPanelControl1_KeyDown);
			this->hkoglPanelControl1->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseDown);
			this->hkoglPanelControl1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseMove);
			this->hkoglPanelControl1->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseWheel);
			// 
			// hScrollBar2
			// 
			this->hScrollBar2->Location = System::Drawing::Point(50, 18);
			this->hScrollBar2->Maximum = 110;
			this->hScrollBar2->Minimum = -100;
			this->hScrollBar2->Name = L"hScrollBar2";
			this->hScrollBar2->Size = System::Drawing::Size(149, 17);
			this->hScrollBar2->TabIndex = 5;
			this->hScrollBar2->Value = 85;
			this->hScrollBar2->Scroll += gcnew System::Windows::Forms::ScrollEventHandler(this, &MyForm::hScrollBar2_Scroll);
			// 
			// Normal_Error
			// 
			this->Normal_Error->AutoSize = true;
			this->Normal_Error->Location = System::Drawing::Point(202, 18);
			this->Normal_Error->Name = L"Normal_Error";
			this->Normal_Error->Size = System::Drawing::Size(26, 12);
			this->Normal_Error->TabIndex = 6;
			this->Normal_Error->Text = L"0.85";
			// 
			// Setting
			// 
			this->Setting->Controls->Add(this->Ignore_Num);
			this->Setting->Controls->Add(this->hScrollBar3);
			this->Setting->Controls->Add(this->face_Ignore);
			this->Setting->Controls->Add(this->Scale_Num);
			this->Setting->Controls->Add(this->ScaleNum);
			this->Setting->Controls->Add(this->label2);
			this->Setting->Controls->Add(this->Face_Error);
			this->Setting->Controls->Add(this->hScrollBar1);
			this->Setting->Controls->Add(this->label3);
			this->Setting->Controls->Add(this->label1);
			this->Setting->Controls->Add(this->hScrollBar2);
			this->Setting->Controls->Add(this->Normal_Error);
			this->Setting->Location = System::Drawing::Point(363, 86);
			this->Setting->Name = L"Setting";
			this->Setting->Size = System::Drawing::Size(238, 154);
			this->Setting->TabIndex = 7;
			this->Setting->TabStop = false;
			this->Setting->Text = L"Setting";
			this->Setting->Visible = false;
			// 
			// Ignore_Num
			// 
			this->Ignore_Num->AutoSize = true;
			this->Ignore_Num->Location = System::Drawing::Point(201, 87);
			this->Ignore_Num->Name = L"Ignore_Num";
			this->Ignore_Num->Size = System::Drawing::Size(11, 12);
			this->Ignore_Num->TabIndex = 16;
			this->Ignore_Num->Text = L"0";
			// 
			// hScrollBar3
			// 
			this->hScrollBar3->Location = System::Drawing::Point(71, 87);
			this->hScrollBar3->Maximum = 1000;
			this->hScrollBar3->Name = L"hScrollBar3";
			this->hScrollBar3->Size = System::Drawing::Size(128, 18);
			this->hScrollBar3->TabIndex = 15;
			this->hScrollBar3->Scroll += gcnew System::Windows::Forms::ScrollEventHandler(this, &MyForm::hScrollBar3_Scroll);
			// 
			// face_Ignore
			// 
			this->face_Ignore->AutoSize = true;
			this->face_Ignore->Location = System::Drawing::Point(9, 93);
			this->face_Ignore->Name = L"face_Ignore";
			this->face_Ignore->Size = System::Drawing::Size(55, 12);
			this->face_Ignore->TabIndex = 14;
			this->face_Ignore->Text = L"faceIgnore";
			// 
			// Scale_Num
			// 
			this->Scale_Num->AutoSize = true;
			this->Scale_Num->Location = System::Drawing::Point(140, 61);
			this->Scale_Num->Margin = System::Windows::Forms::Padding(2, 0, 2, 0);
			this->Scale_Num->Name = L"Scale_Num";
			this->Scale_Num->Size = System::Drawing::Size(11, 12);
			this->Scale_Num->TabIndex = 13;
			this->Scale_Num->Text = L"1";
			// 
			// ScaleNum
			// 
			this->ScaleNum->Location = System::Drawing::Point(64, 59);
			this->ScaleNum->Margin = System::Windows::Forms::Padding(2);
			this->ScaleNum->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 5, 0, 0, 0 });
			this->ScaleNum->Name = L"ScaleNum";
			this->ScaleNum->Size = System::Drawing::Size(58, 22);
			this->ScaleNum->TabIndex = 12;
			this->ScaleNum->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			this->ScaleNum->ValueChanged += gcnew System::EventHandler(this, &MyForm::ScaleNum_ValueChanged);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(7, 61);
			this->label2->Margin = System::Windows::Forms::Padding(2, 0, 2, 0);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(29, 12);
			this->label2->TabIndex = 11;
			this->label2->Text = L"Scale";
			// 
			// Face_Error
			// 
			this->Face_Error->AutoSize = true;
			this->Face_Error->Location = System::Drawing::Point(199, 39);
			this->Face_Error->Name = L"Face_Error";
			this->Face_Error->Size = System::Drawing::Size(26, 12);
			this->Face_Error->TabIndex = 10;
			this->Face_Error->Text = L"0.01";
			// 
			// hScrollBar1
			// 
			this->hScrollBar1->Location = System::Drawing::Point(50, 39);
			this->hScrollBar1->Maximum = 1000;
			this->hScrollBar1->Name = L"hScrollBar1";
			this->hScrollBar1->Size = System::Drawing::Size(149, 17);
			this->hScrollBar1->TabIndex = 9;
			this->hScrollBar1->Value = 1;
			this->hScrollBar1->Scroll += gcnew System::Windows::Forms::ScrollEventHandler(this, &MyForm::hScrollBar1_Scroll);
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(7, 39);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(26, 12);
			this->label3->TabIndex = 8;
			this->label3->Text = L"Face";
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(7, 18);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(40, 12);
			this->label1->TabIndex = 7;
			this->label1->Text = L"Normal";
			// 
			// faceClusterpopulosityToolStripMenuItem
			// 
			this->faceClusterpopulosityToolStripMenuItem->Name = L"faceClusterpopulosityToolStripMenuItem";
			this->faceClusterpopulosityToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->faceClusterpopulosityToolStripMenuItem->Text = L"faceCluster(populosity)";
			this->faceClusterpopulosityToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::faceClusterpopulosityToolStripMenuItem_Click);
			//
			//hypothesismesh
			//
			this->hypothesismeshToolStripMenuItem->Name = L"hypothesismeshToolStripMenuItem";
			this->hypothesismeshToolStripMenuItem->Size = System::Drawing::Size(204, 22);
			this->hypothesismeshToolStripMenuItem->Text = L"Hypothesis Mesh";
			this->hypothesismeshToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::hypothesismeshToolStripMenuItem_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(613, 454);
			this->Controls->Add(this->Setting);
			this->Controls->Add(this->hkoglPanelControl1);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Name = L"MyForm";
			this->Text = L"OpenMesh_EX";
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->Setting->ResumeLayout(false);
			this->Setting->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ScaleNum))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

private: System::Void hkoglPanelControl1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e)
{
	glEnable(GL_COLOR_MATERIAL);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	point center;
	center[0] = 0.0;
	center[1] = 0.0;
	center[2] = 0.0;
	camera.setupGL(xf * center, 1.0);

	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glMultMatrixd((double *)xf);
	if (mesh != NULL)
	{
		if (paintStatus == 0)
			mesh->Render_SolidWireframe();
		else if (paintStatus == 1)
			mesh->Render_newMesh();
	}
	glPopMatrix();
}
private: System::Void hkoglPanelControl1_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left ||
		e->Button == System::Windows::Forms::MouseButtons::Middle)
	{
		point center;
		Mouse_State = Mouse::NONE;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
	}
}
private: System::Void hkoglPanelControl1_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left)
	{
		point center;
		Mouse_State = Mouse::ROTATE;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}

	if (e->Button == System::Windows::Forms::MouseButtons::Middle)
	{
		point center;
		Mouse_State = Mouse::MOVEXY;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
}
private: System::Void hkoglPanelControl1_MouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Delta < 0)
	{
		point center;
		Mouse_State = Mouse::WHEELDOWN;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
	else
	{
		point center;
		Mouse_State = Mouse::WHEELUP;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
}
private: System::Void hkoglPanelControl1_KeyDown(System::Object^ sender, KeyEventArgs^ e) 
{
	if (e->KeyCode == Keys::S) 
	{
		Setting->Visible = !Setting->Visible;
	}
}

private: System::Void loadModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	openModelDialog->Filter = "Model(*.obj)|*obj";
	openModelDialog->Multiselect = false;
	openModelDialog->ShowDialog();
}
private: System::Void openModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
{
	std::string filename;
	MarshalString(openModelDialog->FileName, filename);

	if (mesh != NULL)
		delete mesh;

	mesh = new Tri_Mesh;

	if (ReadFile(filename, mesh))
		std::cout << filename << std::endl;
	mesh->calFaceData();
	mesh->facefaceCount();
	//mesh->contructbbox();
	hkoglPanelControl1->Invalidate();
}
private: System::Void saveModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	saveModelDialog->Filter = "Model(*.obj)|*obj";
	saveModelDialog->ShowDialog();
}
private: System::Void saveModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
{
	std::string filename;
	MarshalString(saveModelDialog->FileName, filename);

	if (SaveFile(filename, mesh->hypmesh))
		std::cout << filename << std::endl;
}

private: System::Void hScrollBar2_Scroll(System::Object^ sender, System::Windows::Forms::ScrollEventArgs^ e) 
{
	Normal_Error->Text = ((float)(hScrollBar2->Value) / 100.0f).ToString();
	setNormalError(((float)(hScrollBar2->Value) / 100.0f));
}
private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) 
{
	mesh->setFaceColor();
}
private: System::Void connectColorToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	if (mesh != NULL)
	{
		std::cout << mesh->facefaceCount();
	}
}
private: System::Void closeFaceToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	closeFace = true;
	//mesh->setFaceColor(closeFace);
	mesh->setCloseFaceColor();
}
private: System::Void nonCloseFaceToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	closeFace = false;
	mesh->setFaceColor();
}
private: System::Void faceClusterToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	mesh->faceCluster(false);
}

private: System::Void faceClusterCloseFaceToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	mesh->faceCluster(true);
}


	private:System::Void faceClusterpopulosityToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
	{
		mesh->facepopCluster();
	}
	private:System::Void hypothesismeshToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
	{
		//mesh->bboxmesh->hypothesismesh();
		 mesh->remesh();
	}

private: System::Void hScrollBar1_Scroll(System::Object^ sender, System::Windows::Forms::ScrollEventArgs^ e) 
{
	Face_Error->Text = ((float)(hScrollBar1->Value) / 100.0f).ToString();
	setFaceError(((float)(hScrollBar1->Value) / 100.0f * getScale()));
}
private: System::Void ScaleNum_ValueChanged(System::Object^ sender, System::EventArgs^ e) 
{
	 Scale_Num->Text = ((float)(pow(10, (int)ScaleNum->Value))).ToString();
	 setScale((float)(pow(10, (int)ScaleNum->Value)));
}
private: System::Void hScrollBar3_Scroll(System::Object^ sender, System::Windows::Forms::ScrollEventArgs^ e) 
{
	Ignore_Num->Text = ((float)(hScrollBar3->Value) / 100.0f * getScale()).ToString();
	setFaceIgnore(((float)(hScrollBar3->Value) / 100.0f * getScale()));
	if (mesh != NULL)mesh->CheckFaceAllow();
}
private: System::Void saveBoundaryToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	if (mesh != NULL && mesh->state == 5) {
		mesh->writeBoundaryFile("boundary");
	}
}
private: System::Void hkoglPanelControl1_Load(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void originToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	paintStatus = 0;
}
private: System::Void newToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	paintStatus = 1;
}
private: System::Void saveNewModelToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	mesh->saveNewModel();
}
private: System::Void saveNewModelMtlToolStripMenuItem_Click(System::Object^ sender, System::EventArgs^ e) 
{
	char* path = "./Assets/cache/";
	char* name = "test";
	mesh->saveNewModelMtl(path, name);
}
};
}
