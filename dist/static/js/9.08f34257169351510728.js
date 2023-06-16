webpackJsonp([9],{1277:function(e,t){e.exports={render:function(){var e=this,t=e.$createElement,a=e._self._c||t;return a("Card",[a("Form",{ref:"formValidate",attrs:{model:e.formValidate,rules:e.ruleValidate,"label-width":80}},[a("Form-item",{attrs:{label:"姓名",prop:"name"}},[a("Input",{attrs:{placeholder:"请输入姓名"},model:{value:e.formValidate.name,callback:function(t){e.$set(e.formValidate,"name ",t)},expression:" formValidate.name "}})],1),e._v(" "),a("Form-item",{attrs:{label:"邮箱",prop:"mail"}},[a("Input",{attrs:{placeholder:"请输入邮箱"},model:{value:e.formValidate.mail,callback:function(t){e.$set(e.formValidate,"mail ",t)},expression:" formValidate.mail "}})],1),e._v(" "),a("Form-item",{attrs:{label:"城市",prop:"city"}},[a("Select",{attrs:{placeholder:"请选择所在地"},model:{value:e.formValidate.city,callback:function(t){e.$set(e.formValidate,"city",t)},expression:"formValidate.city"}},[a("Option",{attrs:{value:"beijing"}},[e._v("北京市")]),e._v(" "),a("Option",{attrs:{value:"shanghai"}},[e._v("上海市")]),e._v(" "),a("Option",{attrs:{value:"shenzhen"}},[e._v("深圳市")])],1)],1),e._v(" "),a("Form-item",{attrs:{label:"生日"}},[a("Row",[a("Col",{attrs:{span:"11"}},[a("Form-item",{attrs:{prop:"date"}},[a("Date-picker",{attrs:{type:"date",placeholder:"选择日期"},model:{value:e.formValidate.date,callback:function(t){e.$set(e.formValidate,"date",t)},expression:"formValidate.date"}})],1)],1)],1)],1),e._v(" "),a("Form-item",{attrs:{label:"性别",prop:"gender"}},[a("Radio-group",{model:{value:e.formValidate.gender,callback:function(t){e.$set(e.formValidate,"gender",t)},expression:"formValidate.gender"}},[a("Radio",{attrs:{label:"male"}},[e._v("男")]),e._v(" "),a("Radio",{attrs:{label:"female"}},[e._v("女")])],1)],1),e._v(" "),a("Form-item",{attrs:{label:"爱好",prop:"interest"}},[a("Checkbox-group",{model:{value:e.formValidate.interest,callback:function(t){e.$set(e.formValidate,"interest",t)},expression:"formValidate.interest"}},[a("Checkbox",{attrs:{label:"吃饭"}}),e._v(" "),a("Checkbox",{attrs:{label:"睡觉"}}),e._v(" "),a("Checkbox",{attrs:{label:"跑步"}}),e._v(" "),a("Checkbox",{attrs:{label:"看电影"}})],1)],1),e._v(" "),a("Form-item",{attrs:{label:"介绍",prop:"desc"}},[a("Input",{attrs:{type:"textarea",autosize:{minRows:2,maxRows:5},placeholder:"请输入..."},model:{value:e.formValidate.desc,callback:function(t){e.$set(e.formValidate,"desc",t)},expression:"formValidate.desc"}})],1),e._v(" "),a("Form-item",[a("Button",{attrs:{type:"primary"},on:{click:function(t){e.handleSubmit("formValidate")}}},[e._v("提交")]),e._v(" "),a("Button",{staticStyle:{"margin-left":"8px"},attrs:{type:"ghost"},on:{click:function(t){e.handleReset("formValidate")}}},[e._v("重置")])],1)],1)],1)},staticRenderFns:[]}},784:function(e,t,a){var r=a(52)(a(966),a(1277),null,null,null);e.exports=r.exports},966:function(e,t,a){"use strict";Object.defineProperty(t,"__esModule",{value:!0}),t.default={data:function(){return{formValidate:{name:"",mail:"",city:"",gender:"",interest:[],date:"",time:"",desc:""},ruleValidate:{name:[{required:!0,message:"姓名不能为空",trigger:"blur"}],mail:[{required:!0,message:"邮箱不能为空",trigger:"blur"},{type:"email",message:"邮箱格式不正确",trigger:"blur"}],city:[{required:!0,message:"请选择城市",trigger:"change"}],gender:[{required:!0,message:"请选择性别",trigger:"change"}],interest:[{required:!0,type:"array",min:1,message:"至少选择一个爱好",trigger:"change"},{type:"array",max:2,message:"最多选择两个爱好",trigger:"change"}],date:[{required:!0,type:"date",message:"请选择日期",trigger:"change"}],desc:[{required:!0,message:"请输入个人介绍",trigger:"blur"},{type:"string",min:20,message:"介绍不能少于20字",trigger:"blur"}]}}},methods:{handleSubmit:function(e){var t=this;this.$refs[e].validate(function(e){e?t.$Message.success("提交成功!"):t.$Message.error("表单验证失败!")})},handleReset:function(e){this.$refs[e].resetFields()}}}}});
//# sourceMappingURL=9.08f34257169351510728.js.map