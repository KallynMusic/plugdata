#pragma once
#include <nanovg.h>
#include "NanoVGGraphicsContext.h"

class NVGComponent
{
    
public:
            
    NVGComponent(Component*);
    
    void renderComponentFromImage(NVGcontext* nvg, Component& component, float scale);
    
    static NVGcolor convertColour(Colour c);
    NVGcolor findNVGColour(int colourId);
    
    virtual void render(NVGcontext*) {};
    
private:
    std::unique_ptr<NanoVGGraphicsContext> nvgLLGC;
    
    struct CachedImage
    {
        int imageId = 0;
        int lastHeight;
        int lastWidth;
    };
    
    CachedImage cachedImage;
    Component& component;
};
