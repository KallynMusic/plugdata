/*
 // Copyright (c) 2021-2022 Timothy Schoen.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#pragma once
#include <utility>
extern "C"
{
#include <pd-lua/lua/lua.h>
#include <pd-lua/lua/lauxlib.h>
#include <pd-lua/lua/lualib.h>
}

#include "Components/BouncingViewport.h"
#include "Object.h"
#include "PluginEditor.h"
#include "Objects/ObjectBase.h"
#include "Sidebar/Sidebar.h"
#include "Components/MarkupDisplay.h"

class CommandProcessor
{
public:
    virtual SmallArray<std::pair<int, String>> executeCommand(pd::Instance* pd, String message) = 0;
};

class LuaExpressionParser {
public:
    using LuaResult = std::variant<double, String>;

    // Constructor: initialize Lua
    LuaExpressionParser(pd::Instance* pd)
        : pd(pd)
    {
        L = luaL_newstate();
        luaL_openlibs(L); // Load Lua standard libraries

        luaopen_math(L);

        // Create the global "pd" table and set up "pd.post"
        lua_newtable(L); // Create a new table for "pd"

        lua_pushlightuserdata(L, this);
        lua_pushcclosure(L, LuaExpressionParser::luaPost, 1);  // 1 upvalue for "pd"
        lua_setfield(L, -2, "post");  // Sets pd.post in the table

        lua_pushlightuserdata(L, this);
        lua_pushcclosure(L, LuaExpressionParser::luaEval, 1); // 1 upvalue for "LuaExpressionParser"
        lua_setfield(L, -2, "eval"); // Sets pd.eval in the table

        lua_setglobal(L, "pd");  // Set the table as the global "pd"
    }

    // Destructor: close Lua
    ~LuaExpressionParser()
    {
        if (L) {
            lua_close(L);
        }
    }

    void executeScript(const juce::String& filePath) {
        // Load the script without executing it
        if (luaL_loadfile(L, filePath.toRawUTF8()) == LUA_OK) {
            // Set the environment of the loaded chunk to _G to make everything global
            lua_pushglobaltable(L);  // Push _G onto the stack
            lua_setupvalue(L, -2, 1);  // Set _G as the environment of the chunk

            // Execute the chunk, which will register all functions globally
            if (lua_pcall(L, 0, 0, 0) != LUA_OK) {
                const char* error = lua_tostring(L, -1);
                pd->logError("Error executing Lua script: " + juce::String::fromUTF8(error));
                lua_pop(L, 1); // Remove error message from stack
            }
        } else {
            const char* error = lua_tostring(L, -1);
            pd->logError("Error loading Lua script: " + juce::String::fromUTF8(error));
            lua_pop(L, 1); // Remove error message from stack
        }
    }


    // Function to execute an expression and return result as LuaResult (either double or string)
    LuaResult executeExpression(String const& expression, bool hasReturnValue)
    {
        String luaCode = "local __eval = function()\n";
        if(hasReturnValue) luaCode += "\nreturn ";
        luaCode += expression.trim(); // Append the expression without altering it
        luaCode += R"(
        end
        local success, result = pcall(__eval)
        if success then
            return result
        else
            error(result)
        end
        )";

        // Run the Lua code and check for errors
        if (luaL_dostring(L, luaCode.toRawUTF8()) == LUA_OK) {
            if (lua_isnumber(L, -1)) {
                double result = lua_tonumber(L, -1);
                lua_pop(L, 1); // Remove result from stack
                return result;
            } else if (lua_isstring(L, -1)) {
                String result = lua_tostring(L, -1);
                lua_pop(L, 1);  // Remove result from stack
                return result;
            } else {
                lua_pop(L, 1); // Remove nil result
                return "";     // Default to empty string if not a number or string
            }
        } else {
            char const* error = lua_tostring(L, -1);
            pd->logError("Lua error: " + String::fromUTF8(error));
            lua_pop(L, 1); // Remove error message from stack
            return "";     // Return empty string on error
        }
    }

    static int luaPost(lua_State* L)
    {
        // Retrieve the LuaWrapper instance (via userdata, upvalues, etc.)
        auto* pd = reinterpret_cast<LuaExpressionParser*>(lua_touserdata(L, lua_upvalueindex(1)))->pd;

        if (lua_isstring(L, 1)) {
            String message = lua_tostring(L, 1);
            pd->logMessage(message); // Use Pd instance to post message
        } else {
            pd->logError("pd.post requires a string argument");
        }
        return 0;
    }
    
    static int luaEval(lua_State* L)
    {
        auto* parser = reinterpret_cast<LuaExpressionParser*>(lua_touserdata(L, lua_upvalueindex(1)));
        
        if (lua_isstring(L, 1)) {
            String command = lua_tostring(L, 1);
            auto result = parser->commandInput->executeCommand(parser->pd, command);  // Execute the command
            
            // Create a Lua table to store result messages
            lua_newtable(L); // Creates an empty table on the stack
            
            int index = 1;
            for (const auto& [type, string] : result) {
                if (type == 0) {
                    lua_pushstring(L, string.toRawUTF8()); // Push the string onto the Lua stack
                    lua_rawseti(L, -2, index);             // Set it in the table at position `index`
                    ++index;
                }
            }
            
            return 1;
            
        }
        
        parser->pd->logError("pd.eval requires a string argument");
        return 0;
    }
    
    void setCommandProcessor(CommandProcessor* newCommandInput)
    {
        commandInput = newCommandInput;
    }

private:
    lua_State* L; // Lua state
    pd::Instance* pd;
    CommandProcessor* commandInput = nullptr;
};

class CommandInput final
    : public Component
    , public KeyListener
    , public CommandProcessor, public MarkupDisplay::URLHandler {
public:
    CommandInput(PluginEditor* editor) : editor(editor)
    {
        // We need to set the target for the command manager, otherwise it will default to PlugDataApp and fail to find CommandID
        editor->commandManager.setFirstCommandTarget(editor);
        // Get the application command id key to toggle show/hide of the command prompt
        // We need to know this as the command prompt gets keyboard focus
        // So the command prompt needs to dismiss itself when the CommandID key is pressed
        auto* keyMappings = editor->commandManager.getKeyMappings();
        auto keyPresses = keyMappings->getKeyPressesAssignedToCommand(CommandIDs::ShowCommandInput);
        commandIDToggleShowKey = keyPresses.getFirst();

        if(!luas.contains(editor->pd))
        {
            luas[editor->pd] = std::make_unique<LuaExpressionParser>(editor->pd);
        }
        lua = luas[editor->pd].get();
        
        updateCommandInputTarget();

        commandInput.setMultiLine(true);
        commandInput.setReturnKeyStartsNewLine(false);
        
        commandInput.onTextChange = [this](){
            currentCommand = commandInput.getText();
            updateSize();
        };

        commandInput.onReturnKey = [this, pd = editor->pd]() {
            currentCommand.clear();

            auto text = commandInput.getText();
            
            if(countBraces(text) > 0)
            {
                commandInput.insertTextAtCaret("\n");
                setConsoleTargetName("lua");
                updateSize();
                return;
            }
            
            auto result = executeCommand(pd, text);
            for(auto& [type, message] : result)
            {
                if(type == 0)
                {
                    pd->logMessage(message);
                }
                else {
                    pd->logError(message);
                }
            }
            auto isUniqueCommand = commandHistory.empty() ? true : commandHistory.front() != commandInput.getText();
            if(!commandInput.isEmpty() && isUniqueCommand) {
                commandHistory.push_front(commandInput.getText());
                currentHistoryIndex = -1;
            }
            commandInput.clear();
            updateCommandInputTarget();
            updateSize();
        };
        
        markupDisplay.setURLHandler(this);
        markupDisplay.setFont(Fonts::getVariableFont());
        markupDisplay.setColour(PlugDataColour::canvasBackgroundColourId, findColour(PlugDataColour::levelMeterBackgroundColourId));
        markupDisplay.setMarkupString(documentationString);
        addChildComponent(&markupDisplay);
        
        helpButton.setWantsKeyboardFocus(false);
        helpButton.setClickingTogglesState(true);
        helpButton.onClick = [this](){
            markupDisplay.setVisible(helpButton.getToggleState());
            updateSize();
        };

        addAndMakeVisible(commandInput);
        addAndMakeVisible(clearButton);
        addAndMakeVisible(helpButton);

        clearButton.setWantsKeyboardFocus(false);
        clearButton.onClick = [this](){
            commandInput.clear();
            currentCommand.clear();
            updateCommandInputTarget();
            updateSize();
        };

        commandInput.setBorder({3, 3, 0, 0});
        commandInput.addKeyListener(this);
        commandInput.addMouseListener(this, false);
        commandInput.setFont(Fonts::getDefaultFont().withHeight(15));

        commandInput.setColour(TextEditor::backgroundColourId, Colours::transparentBlack);
        commandInput.setColour(TextEditor::outlineColourId, Colours::transparentBlack);
        commandInput.setColour(TextEditor::focusedOutlineColourId, Colours::transparentBlack);

        if (currentCommand.isNotEmpty())
            commandInput.setText(currentCommand);
        
        updateSize();
    }
        
    void updateSize()
    {
        auto width = std::clamp(commandInput.getTextWidth() + consoleTargetLength + 30, std::max(250, getWidth()), 400);
        setSize(width, std::max(commandInput.getTextHeight(), 15) + 38 + (helpButton.getToggleState() ? 200 : 0));
    }
        
    int countBraces(String const& text)
    {
        int braceCount = 0;

        for (int i = 0; i < text.length(); ++i)
        {
            juce_wchar currentChar = text[i];

            if (currentChar == '{')
            {
                ++braceCount;  // Increment for each '{'
            }
            else if (currentChar == '}')
            {
                --braceCount;  // Decrement for each '}'
            }
        }

        return braceCount;
    }
    
    void updateCommandInputTarget()
    {
        auto name = String("empty");
        if(auto* cnv = editor->getCurrentCanvas()) {
            auto objects = cnv->getSelectionOfType<Object>();
            if (objects.size() == 1) {
                name = objects[0]->getType();
            } else if (objects.size() > 1){
                name = "(" + String(objects.size()) + " selected)";
            }
        }

        setConsoleTargetName(name);
    }

    static UnorderedMap<String, Object*> getUniqueObjectNames(Canvas* cnv)
    {
        UnorderedMap<String, Object*> result;
        UnorderedMap<String, int> nameCount;
        for (auto* object : cnv->objects) {
            if (!object->gui)
                continue;

            auto type = object->gui->getType();
            if (allAtoms.contains(type)) {
                auto& found = nameCount[type];
                found++;
                result[type + "_" + String(found)] = object;
            } else {
                auto tokens = StringArray::fromTokens(object->gui->getText(), false);
                if (tokens.size()) {
                    auto isGui = allGuis.contains(tokens[0]);
                    auto numArgs = std::min(isGui ? 1 : 2, tokens.size());
                    tokens.removeRange(numArgs, tokens.size() - numArgs);

                    auto uniqueName = tokens.joinIntoString("_");
                    auto& found = nameCount[uniqueName];
                    found++;

                    result[uniqueName + "_" + String(found)] = object;
                }
            }
        }

        return result;
    }
        
    SmallArray<Object*> findObjects(Canvas* cnv, String const& name)
    {
        SmallArray<Object*> found;
        auto names = getUniqueObjectNames(cnv);
        if(name.endsWith("*"))
        {
            auto wildcard = name.upToLastOccurrenceOf("*", false, false);
            for(auto [name, ptr] : names)
            {
                if(name.contains(wildcard))
                {
                    found.add(ptr);
                }
            }
        }
        if (names.contains(name)) {
            found.add(names[name]);
        }
        return found;
    }

    String parseExpressions(String const& message)
    {
        String parsedMessage;
        int startPos = 0;
        
        // if the lua expression is not at the start of the message, we expect a return value
        auto hasReturnValue = !message.startsWith("{");

        while (startPos < message.length()) {
            auto openBrace = message.indexOf(startPos, "{");

            if (openBrace == -1) {
                parsedMessage += message.substring(startPos);
                break;
            }

            parsedMessage += message.substring(startPos, openBrace);

            auto closeBrace = message.indexOf(openBrace, "}");
            if (closeBrace == -1) {
                editor->pd->logError("Unmatched '{' in expression.");
                parsedMessage += message.substring(openBrace); // Append remaining message as-is
                break;
            }

            lua->setCommandProcessor(this);
            String luaExpression = message.substring(openBrace + 1, closeBrace);
            auto result = lua->executeExpression(luaExpression, hasReturnValue);

            if (auto doubleResult = std::get_if<double>(&result)) {
                parsedMessage += String(*doubleResult);
            } else if (auto stringResult = std::get_if<String>(&result)) {
                parsedMessage += *stringResult;
            }

            startPos = closeBrace + 1;
        }

        return parsedMessage;
    }

    SmallArray<std::pair<int, String>> executeCommand(pd::Instance* pd, String message) override
    {
        SmallArray<std::pair<int, String>> result;
        
        message = parseExpressions(message.trim());

        auto tokens = StringArray::fromTokens(message, true);

        // Global or canvas message
        if (!tokens[0].startsWith(";") && (consoleTargetName == ">" || consoleTargetName == "lua >" || tokens[0] == ">" || tokens[0] == "deselect" || tokens[0] == "clear"))
        {
            auto selector = hash(tokens[0]);
            switch (selector) {
            case hash("sel"):
            case hash("select"): {
                if (auto* cnv = editor->getCurrentCanvas()) {
                    if (tokens[1].containsOnly("0123456789")) {
                        int index = tokens[1].getIntValue();
                        if (index >= 0 && index < cnv->objects.size()) {
                            // move the window if it needs to be moved
                            editor->highlightSearchTarget(cnv->objects[index]->getPointer(), true);
                            cnv->updateSidebarSelection();
                        } else {
                            result.add({1, "Object index out of bounds"});
                        }
                    } else {
                        auto objects = findObjects(cnv, tokens[1]);
                        for(auto* object : objects) {
                            cnv->setSelected(object, true);
                            cnv->updateSidebarSelection();
                        }
                        if(objects.empty()) pd->logError("No object found for: " + tokens[1]);
                        // TODO: fix highlighting!
                        //if(objects.size()) editor->highlightSearchTarget(objects[0]->getPointer(), true);
                    }
                }
                updateCommandInputTarget();
                break;
            }
            case hash(">"):
            case hash("deselect"): {
                if (auto* cnv = editor->getCurrentCanvas()) {
                    cnv->deselectAll();
                    cnv->updateSidebarSelection();
                }
                updateCommandInputTarget();
                break;
            }
            case hash("ls"):
            case hash("list"): {
                if (auto* cnv = editor->getCurrentCanvas()) {
                    auto names = getUniqueObjectNames(cnv);
                    for (auto& [name, object] : names) {
                        if (allGuis.contains(object->gui->getType())) {
                            result.add({0, name});
                        } else {
                            result.add({0, name + ": " + object->gui->getText()});
                        }
                    }
                }
                break;
            }
            case hash("find"):
            case hash("search"): {
                if (auto* cnv = editor->getCurrentCanvas()) {
                    auto names = getUniqueObjectNames(cnv);
                    for (auto& [name, object] : names) {
                        auto query = tokens[1];
                        query = query.trimCharactersAtEnd("*"); // No need for wildcards here
                        auto text = object->gui->getText();
                        if (text.contains(query) || name.contains(query)) {
                            if (allGuis.contains(object->gui->getType())) {
                                result.add({0, name});
                            } else {
                                result.add({0, name + ": " + object->gui->getText()});
                            }
                        }
                    }
                }

                break;
            }
            case hash("reset"): {
                // Reset lua context
                luas[editor->pd] = std::make_unique<LuaExpressionParser>(editor->pd);
                lua = luas[editor->pd].get();
                break;
            }
            case hash("clear"): {
                commandHistory.clear();
                editor->sidebar->clearConsole();
                if (auto* cnv = editor->getCurrentCanvas()) {
                    cnv->deselectAll();
                    cnv->updateSidebarSelection();
                }
                updateCommandInputTarget();
                break;
            }
            case hash("cnv"):
            case hash("canvas"):
            {
                if(auto* cnv = editor->getCurrentCanvas())
                {
                    auto patchPtr = cnv->patch.getPointer();
                    if (patchPtr && tokens.size() == 1 && tokens[1].containsOnly("0123456789-e.")) {
                        pd->sendDirectMessage(patchPtr.get(), tokens[1].getFloatValue());
                    } else if (patchPtr && tokens.size() == 1) {
                        pd->sendDirectMessage(patchPtr.get(), tokens[1], {});
                    } else if(patchPtr) {
                        SmallArray<pd::Atom> atoms;
                        for (int i = 2; i < tokens.size(); i++) {
                            if (tokens[i].containsOnly("0123456789-e.")) {
                                atoms.add(pd::Atom(tokens[i].getFloatValue()));
                            } else {
                                atoms.add(pd::Atom(pd->generateSymbol(tokens[i])));
                            }
                        }
                        pd->sendDirectMessage(patchPtr.get(), tokens[1], std::move(atoms));
                    }
                    cnv->patch.deselectAll();
                }
                else {
                    pd->logError("No canvas open");
                }
            }
            case hash("script"):
            {
                auto script = pd::Library::findFile(tokens[1] + ".lua");
                if(script.existsAsFile()) {
                    lua->executeScript(script.getFullPathName());
                }
                else {
                    pd->logError("Script not found");
                }
                break;
            }
            case hash("man"):
            {
                switch(hash(tokens[1]))
                {
                    case hash("man"):
                        pd->logMessage("Prints manual for command. Usage: man <command>");
                        break;
                        
                    case hash("?"):
                    case hash("help"):
                        pd->logMessage(tokens[2] + ": Show help");
                        break;
                        
                    case hash("script"):
                        pd->logMessage(tokens[2] + ": Excute a Lua script from your search path. Usage: script <filename>");
                        break;
                        
                    case hash("cnv"):
                    case hash("canvas"):
                        pd->logMessage(tokens[2] + ": Send a message to current canvas. Usage: " + tokens[2] + " <message>");
                        break;
                        
                    case hash("clear"):
                        pd->logMessage(tokens[2] + ": Clear console and command history");
                        break;
                        
                    case hash("reset"):
                        pd->logMessage(tokens[2] + ": Reset Lua interpreter state");
                        break;
                        
                    case hash("sel"):
                    case hash("select"):
                        pd->logMessage(tokens[2] + ": Select an object by ID or index. After selecting objects, you can send messages to them. Usage: " + tokens[2] + " <id> or " + tokens[2] + " <index>");
                        break;
                        
                    case hash(">"):
                    case hash("deselect"):
                        pd->logMessage(tokens[2] + ": Deselects all on current canvas");
                        break;
                        
                    case hash("ls"):
                    case hash("list"):
                        pd->logMessage(tokens[2] + ": Print a list of all object IDs on current canvas");
                        break;
                        
                    case hash("find"):
                    case hash("search"):
                        pd->logMessage(tokens[2] + ": Search object IDs on current canvas. Usage: " + tokens[2] + " <id>.");
                        break;
                }
            }
            case hash("?"):
            case hash("help"):
            {
                helpButton.setToggleState(true, sendNotification);
                break;
            }
            default: {
                // Match a  "name > message" pattern
                if(tokens.size() >= 2 && tokens[1] == ">")
                {
                    auto target = tokens[0];
                    if(tokens.size() == 2)
                    {
                        if(auto* cnv = editor->getCurrentCanvas()) {
                            auto objects = findObjects(cnv, target);
                            for(auto* object : objects) {
                                cnv->setSelected(object, true);
                                cnv->updateSidebarSelection();
                            }
                            if(objects.empty()) pd->logError("No object found for: " + tokens[1]);
                        }
                        break;
                    }
                    
                    tokens.removeRange(0, 2);

                    if(auto* cnv = editor->getCurrentCanvas()) {
                        auto objects = findObjects(cnv, target);
                        for(auto* object : objects) {
                            if(auto* cnv = editor->getCurrentCanvas())
                            {
                                auto objPtr = object->getPointer();
                                if (objPtr && tokens.size() == 1 && tokens[0].containsOnly("0123456789-e.")) {
                                    pd->sendDirectMessage(objPtr, tokens[0].getFloatValue());
                                } else if (objPtr && tokens.size() == 1) {
                                    pd->sendDirectMessage(objPtr, tokens[0], {});
                                } else if(objPtr) {
                                    SmallArray<pd::Atom> atoms;
                                    for (int i = 1; i < tokens.size(); i++) {
                                        if (tokens[i].containsOnly("0123456789-e.")) {
                                            atoms.add(pd::Atom(tokens[i].getFloatValue()));
                                        } else {
                                            atoms.add(pd::Atom(pd->generateSymbol(tokens[i])));
                                        }
                                    }
                                    pd->sendDirectMessage(objPtr, tokens[0], std::move(atoms));
                                }
                            }
                        }
                        if(objects.empty()) pd->logError("No object found for: " + tokens[1]);
                    }
                }
                
                if(!tokens.size()) break;
                tokens.getReference(0) = tokens[0].trimCharactersAtStart(";");
                SmallArray<pd::Atom> atoms;
                for (int i = 2; i < tokens.size(); i++) {
                    if (tokens[i].containsOnly("0123456789-e.")) {
                        atoms.add(pd::Atom(tokens[i].getFloatValue()));
                    } else {
                        atoms.add(pd::Atom(pd->generateSymbol(tokens[i])));
                    }
                }
                pd->sendMessage(tokens[0].toRawUTF8(), tokens[1].toRawUTF8(), std::move(atoms));
                break;
            }
            }
        } else { // object message
            if (auto* cnv = editor->getCurrentCanvas()) {
                for (auto* obj : cnv->getSelectionOfType<Object>()) {
                    if (tokens.size() == 1 && tokens[0].containsOnly("0123456789-e.")) {
                        pd->sendDirectMessage(obj->getPointer(), tokens[0].getFloatValue());
                    } else if (tokens.size() == 1) {
                        pd->sendDirectMessage(obj->getPointer(), tokens[0], {});
                    } else {
                        SmallArray<pd::Atom> atoms;
                        for (int i = 1; i < tokens.size(); i++) {
                            if (tokens[i].containsOnly("0123456789-e.")) {
                                atoms.add(pd::Atom(tokens[i].getFloatValue()));
                            } else {
                                atoms.add(pd::Atom(pd->generateSymbol(tokens[i])));
                            }
                        }
                        pd->sendDirectMessage(obj->getPointer(), tokens[0], std::move(atoms));
                    }
                    
                }
            }
        }
                    
       return result;
    }

    ~CommandInput() override
    {
        onDismiss();
    };
    
    void handleURL(String const& url) override // when documentation links or codeblocks are clicked
    {
        commandInput.setText(url);
    }

    void paintOverChildren(Graphics& g) override
    {
        auto bounds = getLocalBounds().withTrimmedTop(26);
        
        if(helpButton.getToggleState())
        {
            bounds.removeFromTop(200);
        }
        
        g.setColour(findColour(PlugDataColour::dataColourId));
        g.setFont(Fonts::getSemiBoldFont().withHeight(15));
        g.drawText(consoleTargetName, bounds.getX() + 7, bounds.getY(), consoleTargetLength, bounds.getHeight() - 3, Justification::centredLeft);
    }
        
    void paint(Graphics& g) override
    {
        auto bounds = getLocalBounds();
        g.setFont(Fonts::getSemiBoldFont().withHeight(15));
        g.setColour(findColour(PlugDataColour::panelTextColourId));
        g.drawText("Command input", bounds.removeFromTop(22), Justification::centred);
        
        bounds.removeFromTop(4);
        
        if(helpButton.getToggleState())
        {
            bounds.removeFromTop(200);
        }
        
        g.setColour(findColour(PlugDataColour::levelMeterBackgroundColourId));
        g.fillRoundedRectangle(bounds.reduced(2, 2).toFloat(), Corners::defaultCornerRadius);
    }

    void resized() override
    {
        auto inputBounds = getLocalBounds().withTrimmedTop(26);
        if(helpButton.getToggleState())
        {
            markupDisplay.setBounds(inputBounds.removeFromTop(196));
            inputBounds.removeFromTop(4);
        }
        commandInput.setBounds(inputBounds.withTrimmedLeft(consoleTargetLength).withTrimmedRight(30));
        auto buttonBounds = inputBounds.removeFromRight(30);
        clearButton.setBounds(buttonBounds);
        helpButton.setBounds(getLocalBounds().removeFromTop(22).removeFromRight(22));
    }

    void setConsoleTargetName(String const& target)
    {
        consoleTargetName = target + " >";
        if (target == "empty")
            consoleTargetName = ">";
        consoleTargetLength = CachedStringWidth<15>::calculateStringWidth(consoleTargetName) + 4;
        commandInput.setBounds(commandInput.getBounds().withLeft(consoleTargetLength));
        repaint();
    }

    void setHistoryCommand()
    {
        if (currentHistoryIndex < 0) {
            commandInput.setText("");
            currentHistoryIndex = -1;
        } else if (currentHistoryIndex < commandHistory.size()) {
            auto command = commandHistory[currentHistoryIndex];
            auto isMultiLine = command.containsChar('\n');
            if(isMultiLine) setConsoleTargetName("lua");
            else updateCommandInputTarget();
            commandInput.setText(command);
        } else {
            currentHistoryIndex = commandHistory.size() - 1;
        }
    }

    bool keyPressed(KeyPress const& key, Component*) override
    {
        if (key.getKeyCode() == KeyPress::returnKey && key.getModifiers().isShiftDown()) {
            commandInput.insertTextAtCaret("\n");
            updateSize();
            return true;
        }
        else if (key.getKeyCode() == KeyPress::upKey && !commandInput.getText().containsChar('\n')) {
            currentHistoryIndex++;
            setHistoryCommand();
            return true;
        }
        else if (key.getKeyCode() == KeyPress::downKey && !commandInput.getText().containsChar('\n')) {
            currentHistoryIndex--;
            setHistoryCommand();
            return true;
        }
        else if (key.getKeyCode() == KeyPress::escapeKey) {
            if (auto* cnv = editor->getCurrentCanvas())
                cnv->deselectAll();
            updateCommandInputTarget();
            return true;
        }
        else if (key.getKeyCode() == commandIDToggleShowKey.getKeyCode()) {
            editor->commandManager.invokeDirectly(CommandIDs::ShowCommandInput, false);
           return true;
        }
        return false;
    }

    static std::deque<String>& getCommandHistory()
    {
        return commandHistory;
    }

    static void setCommandHistory(StringArray& commands)
    {
        commandHistory.clear();

        for (auto& command : commands){
            commandHistory.push_back(command);
        }
    };

    private:
    PluginEditor* editor;
    static inline UnorderedMap<pd::Instance*, std::unique_ptr<LuaExpressionParser>> luas = UnorderedMap<pd::Instance*, std::unique_ptr<LuaExpressionParser>>();
    LuaExpressionParser* lua;
    
    int consoleTargetLength = 10;
    String consoleTargetName = ">";

    int currentHistoryIndex = -1;
    static inline std::deque<String> commandHistory;
    static inline String currentCommand;

    KeyPress commandIDToggleShowKey;

    TextEditor commandInput;
    SmallIconButton clearButton = SmallIconButton(Icons::ClearText);
    SmallIconButton helpButton = SmallIconButton(Icons::Help);
    
    MarkupDisplay::MarkupDisplayComponent markupDisplay;
        
    static inline String documentationString = {
        "Command input allows you to quickly send commands to objects, pd or the canvas.\n"
        "The following commands are available:\n"
        "\n"
        "- man <command>: print manual for command\n"
        "- list/ls: list all object IDs in the current canvas\n"
        "- search: search for an object ID in current canvas\n"
        "- select/sel <id>: selects an object\n"
        "- deselect: deselect all objects\n"
        "- clear: clears console and command state\n"
        "- reset: clear lua state\n"
        "- canvas: send message to canvas\n"
        "    - canvas obj <x> <y> <name>: create text object\n"
        "    - canvas msg <x> <y> <name>: create message object\n"
        "- script: load lua script from search path\n"
        "- pd: send a message to pd, for example:\n"
        "    - pd dsp <int>: set DSP state\n"
        "    - pd pluginmode: enter plugin mode\n"
        "    - pd quit: quit plugdata\n"
        "- <id> > <message>: sends message to object, for example:\n"
        "    ```tgl_1 > 1```\n to send a 1 to toggle object 1\n"
        "\n"
        "Once you have selected an object, all messages you send become direct messages to that object. So to send a float to a toggle, select it, and enter \"1\". You can deselect an object with \"deselect\", or the shorthand \">\"\n"
        "\n"
        "You can also do dynamic patching with the \"canvas\" command.\n"
        "\n"
        "```canvas obj 20 50 metro 200```\n"
        "\n"
        "This will create a \"metro 200\" object at the 20,50 coordinate.\n"
        "\n"
        "\n"
        "You can also use Lua expressions to generate values, or automate tasks. Lua expressions are written inside brackets. For example, to randomise the colour of a toggle, select a toggle and enter:\n"
        "\n"
        "```color {math.random() * 200}```\n"
        "\n"
        "Lua can also call back to the command input by calling pd.eval(\"your command here\"). For example:\n"
        "\n"
        "```{ pd.eval(\"sel tgl_1\") }```\n"
        "\n"
        "Will select the first toggle from within lua, which you could then send more messages to from Lua\n"
        "\n"
        "To write multi-line lua expressions, leave an open bracket and hit enter, so you could write:\n"
        "\n"
        "```{\n"
        "for i = 1, 20 do\n"
        "    for j = 1, 20 do\n"
        "        pd.eval(\"canvas obj \" .. tostring(i * 28) .. \" \" .. tostring(j * 28) .. \" tgl\")\n"
        "    end\n"
        "end\n"
        "}```\n"
        "\n"
        "To generate a 20x20 grid of toggle objects."
    };

public:

    std::function<void()> onDismiss = [](){};
        
    static inline const UnorderedSet<String> allAtoms = { "floatbox", "symbolbox", "listbox", "gatom" };

    static inline const UnorderedSet<String> allGuis = { "bng",
                                                   "hsl",
                                                   "vsl",
                                                   "slider",
                                                   "tgl",
                                                   "nbx",
                                                   "vradio",
                                                   "hradio",
                                                   "vu",
                                                   "cnv",
                                                   "keyboard",
                                                   "pic",
                                                   "scope~",
                                                   "function",
                                                   "note",
                                                   "knob",
                                                   "message",
                                                   "comment",
                                                   "canvas",
                                                   "bicoeff",
                                                   "messbox",
                                                   "pad",
                                                   "button"};
};
