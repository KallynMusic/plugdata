/*
 // Copyright (c) 2024 Timothy Schoen
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

class BorderResizer : public Component
    , public NVGComponent
    , public Value::Listener {
    ComponentDragger dragger;
    Canvas* cnv;

public:
    std::function<void()> onDrag = []() { };

    BorderResizer(Canvas* parentCanvas)
        : NVGComponent(this)
        , cnv(parentCanvas)
    {
        setSize(8, 8);
        setRepaintsOnMouseActivity(true);

        cnv->patchHeight.addListener(this);
        cnv->patchWidth.addListener(this);
    }

    ~BorderResizer()
    {
        cnv->patchHeight.removeListener(this);
        cnv->patchWidth.removeListener(this);
    }

    void mouseDown(MouseEvent const& e) override
    {
        if (cnv->showBorder) {
            dragger.startDraggingComponent(this, e);
        }
    }

    void mouseDrag(MouseEvent const& e) override
    {
        if (cnv->showBorder) {
            auto constrainedPoint = getLocalPoint(cnv, Rectangle<int>(cnv->canvasOrigin.x + 11, cnv->canvasOrigin.y + 11, cnv->canvasOrigin.x, cnv->canvasOrigin.y).getConstrainedPoint(e.getEventRelativeTo(cnv).getPosition()));
            dragger.dragComponent(this, e.withNewPosition(constrainedPoint), nullptr);

            cnv->patchHeight.removeListener(this);
            cnv->patchWidth.removeListener(this);
            onDrag();
        }
    }

    void mouseUp(MouseEvent const& e) override
    {
        cnv->patchHeight.addListener(this);
        cnv->patchWidth.addListener(this);
    }

    void render(NVGcontext* nvg) override
    {
        NVGScopedState state(nvg);
        nvgSave(nvg);
        nvgTranslate(nvg, getX(), getY());
        auto bounds = getLocalBounds().reduced(isMouseOver() ? 0 : 2);
        nvgDrawRoundedRect(nvg, bounds.getX(), bounds.getY(), bounds.getWidth(), bounds.getHeight(), findNVGColour(PlugDataColour::canvasDotsColourId), findNVGColour(PlugDataColour::canvasDotsColourId), bounds.getWidth() / 2.0f);
        nvgRestore(nvg);
    }

    void valueChanged(Value& v) override
    {
        if (v.refersToSameSourceAs(cnv->patchHeight)) {
            setCentrePosition(getBoundsInParent().getCentreX(), getValue<int>(cnv->patchHeight) + cnv->canvasOrigin.y);
        } else if (v.refersToSameSourceAs(cnv->patchWidth)) {
            setCentrePosition(getValue<int>(cnv->patchWidth) + cnv->canvasOrigin.x, getBoundsInParent().getCentreY());
        }
    }
};
