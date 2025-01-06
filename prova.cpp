#include <TCanvas.h>
#include <TH1D.h>
#include <TApplication.h>

int main() {
    TApplication theApp("App", 0, 0); // Crea un'applicazione ROOT
    
    // Crea un canvas
    TCanvas* canvas = new TCanvas("c1", "Canvas per Istogramma", 800, 600);
    
    // Crea e disegna un istogramma
    TH1D* h = new TH1D("h", "Istogramma di esempio", 100, 0, 10);
    h->FillRandom("gaus", 1000);
    h->Draw();
    
    // Mostra il canvas
    canvas->Update();
    
    // Esegui l'applicazione ROOT
    theApp.Run();
    
    return 0;
}
